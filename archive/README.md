# Code and results for the article

[Download the code and results](./link-protection-code-and-results.tar.gz) for *The value of topology information for link protection in random recursive trees*. The 17 files in the download match the [licensed source revision](https://github.com/mgalazka84/recursive-tree-link-protection/commit/97ca5258ab43d609915c24d5df3f5798f4fa729a), excluding the archive directory itself. The C++ and Python files in `code/` and the root `Makefile` are covered by the [MIT License](../LICENSE). A copy of the licence is included in the download.

The source code, build script and numerical results, including the inputs for Table 4, are byte-for-byte identical to those in the [computational revision cited in the article](https://github.com/mgalazka84/recursive-tree-link-protection/commit/66903e22b4858570712db887b280353501a3b996).

The SHA-256 checksum of the compressed file is `4b7aff6d9f63212f78021920dee02a2c0a7688acd0daf698752116eea0a4b1b7`. From the directory containing the download, run the following commands with GNU Make, a C++17 compiler and Python 3 with matplotlib.

~~~sh
sha256sum link-protection-code-and-results.tar.gz
tar -xzf link-protection-code-and-results.tar.gz
cd link-protection
make verify
make reproduce
~~~

Compare the output of the first command with the checksum above. The last two commands check the exact small-tree cases, write the numerical results to `results/`, and regenerate the three figures in `figures/`.
