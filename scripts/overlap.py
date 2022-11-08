#!/usr/bin/env python3
import sys
from difflib import SequenceMatcher

def overlap(s1, s2):
    s = SequenceMatcher(None, s1, s2)
    pos1, _, size = s.find_longest_match()
    return s1[pos1:pos1+size]

if __name__ == "__main__":
    print(overlap(sys.argv[1], sys.argv[2]))
