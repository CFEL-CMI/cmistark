echo "\n*** Line count ***"
find lib/*.* scripts/*.* user-install.csh setup.py MANIFEST [A-LNZ]* -type f | xargs wc -l
echo "\n*** Characters in archive ***"
wc -c dist/cmistark-*.tar.gz
echo "\n"
