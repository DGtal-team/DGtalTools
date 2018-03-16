#!/bin/bash

set -ev
pwd
echo $HOME
rsync -azv --delete --delete-after -e 'ssh -oStrictHostKeyChecking=no -i  /home/travis/build/DGtal-team/DGtalTools/.travis/dgtal_rsa' html/ dgtal@connect.liris.cnrs.fr:/home-projets/dgtal/public_html/doc/tools/nightly


