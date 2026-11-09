# Archived computational materials

The archive contains all 15 tracked files (source, input-free experiments, numerical results and reproduction instructions) from the [computational revision](https://github.com/mgalazka84/recursive-tree-link-protection/commit/aacde39cae7baf8da6114dcbd05010194b6cc204). Its directory has the intrinsic Software Hash Identifier `swh:1:dir:3f82abf5f50e83e629a0f16738f9eeec84085494`. Its compressed bytes have SHA-256 `41568a5b82637d956d2a800970a7619ace3e1a1ddd25e70892c19357906e9ed8`.

To verify the directory and reproduce the calculations, download the archive and run:

~~~sh
tar -xzf link-protection-code-and-results.tar.gz
cd link-protection
git init -q
git add -f -- .
git write-tree
make verify
make reproduce
~~~

The output of `git write-tree` should be `3f82abf5f50e83e629a0f16738f9eeec84085494`. Compute it before running the build, in a freshly extracted directory. The [SWHID specification](https://www.swhid.org/) describes the persistent, content-derived directory identifier.
