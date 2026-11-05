#!/usr/bin/env python3
"""Prepare isolated compile-time BatMIS checkpoint-spacing experiments.

All generated sources, binaries, and indexes go under --build-root. 512/1024
variants compose the unchanged BWTDecode primitive over aligned <=128-row
segments and use a direct packed-BWT checkpoint builder because the historical
packed-byte builder only supports <=256-row minor intervals.
"""
import argparse
import os
import pathlib
import shutil
import subprocess
import time


SPACINGS = (64, 128, 256, 512, 1024)
CFLAGS = "-w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -g -O2"
C_SOURCES = ("BWT.c", "BWTConstruct.c", "MiscUtilities.c", "MemManager.c",
             "TextConverter.c", "QSufSort.c", "DNACount.c")


def run(command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def prepare_sources(repo, root, spacing):
    d = root / str(spacing)
    for sub in ("index", "src", "obj"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    for suffix in ("*.c", "*.h", "*.cpp"):
        for p in (repo / "index").glob(suffix):
            shutil.copy2(p, d / "index" / p.name)
    for name in ("exact_search.cpp", "exact_search.hpp", "batman_zero_structured.cpp"):
        shutil.copy2(repo / "src" / name, d / "src" / name)

    header = d / "index/BWT.h"
    text = header.read_text()
    text, n = __import__("re").subn(r"#define OCC_INTERVAL\s+256",
                                    f"#define OCC_INTERVAL\t\t\t\t{spacing}", text, count=1)
    if n != 1:
        raise RuntimeError(f"Could not set OCC_INTERVAL={spacing}")
    header.write_text(text)

    if spacing > 256:
        bwt = d / "index/BWT.c"
        text = bwt.read_text()
        old = """if (occIndex != index) {
decodeValue = BWTDecode(bwt, occIndex, index, character);
r = -(occIndex > index);
return occValue + (decodeValue & ~r) - (decodeValue & r);
} else {
return occValue;
}"""
        new = """if (occIndex != index) {
unsigned int lowIndex = occIndex < index ? occIndex : index;
unsigned int highIndex = occIndex < index ? index : occIndex;
unsigned int cursor = lowIndex;
decodeValue = 0;
while (cursor < highIndex) {
    unsigned int next;
    if (cursor % CHAR_PER_128 != 0) {
        next = (cursor + CHAR_PER_128 - 1) / CHAR_PER_128 * CHAR_PER_128;
    } else {
        next = cursor + CHAR_PER_128 * 2;
    }
    if (next > highIndex) next = highIndex;
    decodeValue += BWTDecode(bwt, cursor, next, character);
    cursor = next;
}
r = -(occIndex > index);
return occValue + (decodeValue & ~r) - (decodeValue & r);
} else {
return occValue;
}"""
        if old not in text:
            raise RuntimeError("Could not find BWTOccValue decode block")
        bwt.write_text(text.replace(old, new, 1))

        # Full per-checkpoint counts avoid the legacy builder's byte-lane
        # overflow special case, which only adds exactly 256 occurrences.
        construct = d / "index/BWTConstruct.c"
        text = construct.read_text()
        start = text.index("void BWTGenerateOccValueFromBwt(")
        end = text.index("void BWTSaveBwtCodeAndOcc(", start)
        generic_builder = """void BWTGenerateOccValueFromBwt(const unsigned int* bwt, unsigned int* __restrict occValue, unsigned int* __restrict occValueMajor, const unsigned int textLength, const unsigned int* decodeTable) {
unsigned int numberOfValues = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
unsigned int intervalsPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
unsigned int counts[ALPHABET_SIZE] = {0, 0, 0, 0};
unsigned int majorCounts[ALPHABET_SIZE] = {0, 0, 0, 0};
unsigned int pos = 0, checkpoint, c, base, minor, word;
(void)decodeTable;
memset(occValue, 0, BWTOccValueMinorSizeInWord(textLength) * sizeof(unsigned int));
memset(occValueMajor, 0, BWTOccValueMajorSizeInWord(textLength) * sizeof(unsigned int));
for (checkpoint = 0; checkpoint < numberOfValues; ++checkpoint) {
    unsigned int target = checkpoint * OCC_INTERVAL;
    while (pos < target) {
        c = (bwt[pos / CHAR_PER_WORD] >> (BITS_IN_WORD - BIT_PER_CHAR * (pos % CHAR_PER_WORD + 1))) & 3;
        ++counts[c];
        ++pos;
    }
    if (checkpoint % intervalsPerMajor == 0) {
        for (base = 0; base < ALPHABET_SIZE; ++base) {
            majorCounts[base] = counts[base];
            occValueMajor[(checkpoint / intervalsPerMajor) * ALPHABET_SIZE + base] = counts[base];
        }
    }
    for (base = 0; base < ALPHABET_SIZE; ++base) {
        minor = counts[base] - majorCounts[base];
        word = (checkpoint / OCC_VALUE_PER_WORD) * ALPHABET_SIZE + base;
        if (checkpoint % OCC_VALUE_PER_WORD == 0) occValue[word] = minor << BITS_PER_OCC_VALUE;
        else occValue[word] |= minor;
    }
}
}
"""
        construct.write_text(text[:start] + generic_builder + text[end:])
    return d


def compile_variant(d, spacing, repo):
    obj = d / "obj"
    for source in C_SOURCES:
        run(["gcc", *CFLAGS.split(), "-I", str(d / "index"), "-c",
             str(d / "index" / source), "-o", str(obj / (pathlib.Path(source).stem + ".o"))])
    archive = obj / "libfmindex_core.a"
    if archive.exists(): archive.unlink()
    run(["ar", "rcs", str(archive), *[str(obj / (pathlib.Path(s).stem + ".o")) for s in C_SOURCES]])
    run(["g++", "-std=c++17", "-fpermissive", *CFLAGS.split(), "-I", str(d / "index"),
         str(d / "src/batman_zero_structured.cpp"), str(d / "src/exact_search.cpp"),
         str(archive), "-lm", "-o", str(d / "batman_zero_structured")])
    run(["g++", "-std=c++17", *CFLAGS.split(), "-I", str(d / "index"),
         str(d / "index/bwt_build_index.cpp"), str(archive), "-lm", "-o", str(d / "bwt_build_index")])
    for name, binary in (("check_occ_spacing_rank.cpp", "check_occ_spacing_rank"),
                         ("replay_occ_spacing_rank.cpp", "rank_replay")):
        run(["g++", "-std=c++17", "-O2", "-funroll-loops", "-maccumulate-outgoing-args", "-msse2",
             "-I", str(d / "index"), str(repo / "benchmark" / name), str(archive), "-lm",
             "-o", str(d / binary)])


def trace_source(d):
    p = d / "src/exact_search.cpp"
    text = p.read_text().replace("#include <cassert>", "#include <cassert>\n#include <cstdlib>\n#include <cstdint>", 1)
    old = """constexpr int kBranchThreshold = 0;
}"""
    new = """constexpr int kBranchThreshold = 0;
struct RankTraceEntry { std::uint32_t row; std::uint8_t base; std::uint8_t origin; std::uint16_t reserved; };
FILE* gRankTrace = nullptr;
unsigned int tracedRank(const BWT* index, unsigned long row, unsigned int base, std::uint8_t origin) {
    if (gRankTrace) {
        RankTraceEntry entry{static_cast<std::uint32_t>(row), static_cast<std::uint8_t>(base), origin, 0};
        if (std::fwrite(&entry, sizeof(entry), 1, gRankTrace) != 1) std::abort();
    }
    return BWTOccValue(index, static_cast<unsigned int>(row), base);
}
}"""
    if old not in text: raise RuntimeError("Could not add rank trace helper")
    text = text.replace(old, new, 1)
    old = """    lookupSize_ = readLength < 28 ? 3 : 6;

    const int seedCount"""
    new = """    lookupSize_ = readLength < 28 ? 3 : 6;
    if (!gRankTrace) {
        const char* tracePath = std::getenv("BATMIS_RANK_TRACE_PATH");
        if (tracePath) {
            gRankTrace = std::fopen(tracePath, "wb");
            if (!gRankTrace) return false;
        }
    }

    const int seedCount"""
    if old not in text: raise RuntimeError("Could not add trace environment handling")
    text = text.replace(old, new, 1)
    start = text.index("bool ExactSearch::searchExact(")
    end = text.index("bool ExactSearch::search(", start)
    body = text[start:end]
    body = body.replace("BWTOccValue(index_, range.Start, now)", "tracedRank(index_, range.Start, now, 1)")
    body = body.replace("BWTOccValue(index_, range.Start, base)", "tracedRank(index_, range.Start, base, 2)")
    body = body.replace("BWTOccValue(index_, range.End + 1, base)", "tracedRank(index_, range.End + 1, base, 2)")
    text = text[:start] + body + text[end:]
    p.write_text(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--build-root", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--reads", required=True)
    a = ap.parse_args()
    repo = pathlib.Path(a.repo).resolve()
    root = pathlib.Path(a.build_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    dirs = {}
    for spacing in SPACINGS:
        d = prepare_sources(repo, root, spacing)
        compile_variant(d, spacing, repo)
        dirs[spacing] = d

    for spacing, d in dirs.items():
        start = time.monotonic()
        run([str(d / "bwt_build_index"), "-o", str(d), a.reference],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"built_occ_spacing={spacing} index_build_wall_s={time.monotonic()-start:.6f}")
        run([str(d / "check_occ_spacing_rank"), str(d / "revref.fa.bwt"), str(d / "revref.fa.fmv")])

    # The trace is generated once with the 256-row search traversal, before
    # any variant timing; it is then replayed unchanged for every spacing.
    base = dirs[256]
    trace_source(base)
    run(["g++", "-std=c++17", "-fpermissive", *CFLAGS.split(), "-I", str(base / "index"),
         str(base / "src/batman_zero_structured.cpp"), str(base / "src/exact_search.cpp"),
         str(base / "obj/libfmindex_core.a"), "-lm", "-o", str(base / "trace_generator")])
    trace = root / "rank_trace.bin"
    env = os.environ.copy()
    env["BATMIS_RANK_TRACE_PATH"] = str(trace)
    run([str(base / "trace_generator"), "-q", a.reads, "-g", str(base / "ref.fa"),
         "-m", "0", "--quiet", "--timing"], env=env, stdout=subprocess.DEVNULL)
    print(f"rank_trace={trace} bytes={trace.stat().st_size}")


if __name__ == "__main__":
    main()
