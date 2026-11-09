# Code and results for the article

[Download the code and results](./link-protection-code-and-results.tar.gz) for *The value of topology information for link protection in random recursive trees*. The bundle contains the 16 computational files from [revision 66903e22b4858570712db887b280353501a3b996](https://github.com/mgalazka84/recursive-tree-link-protection/commit/66903e22b4858570712db887b280353501a3b996), including the inputs and uncertainty calculations for Table 4. The archive directory itself is excluded.

The SHA-256 checksum of the compressed file is `8f794a17169451a8bf687615118b1d12413ab979df3ad682046dc0e6fa6d4c85`. From the directory containing the download, run the following commands with GNU Make, a C++17 compiler and Python 3 with matplotlib.

~~~sh
sha256sum link-protection-code-and-results.tar.gz
tar -xzf link-protection-code-and-results.tar.gz
cd link-protection
make verify
make reproduce
~~~

Compare the output of the first command with the checksum above. The last two commands check the exact small-tree cases, write the numerical results to `results/`, and regenerate the three figures in `figures/`.
