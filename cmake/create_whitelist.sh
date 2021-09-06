#!/bin/zsh
cd $1;
foreach i (`git diff --name-only HEAD HEAD~5 | grep cpp`)
 echo -e "$WL$i:r:t\c"
 WL=";"
end

