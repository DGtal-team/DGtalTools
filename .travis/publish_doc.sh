#!/bin/bash

set -ev

cd build/
echo $HOME
rsync -azv --delete --delete-after -e 'ssh -oStrictHostKeyChecking=no -i  /home/travis/build/DGtal-team/DGtal/.travis/dgtal_rsa' html/ dgtal@liris.cnrs.fr:/home/dgtal/public_html/doc/tools/


