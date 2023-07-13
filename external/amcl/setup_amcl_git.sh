#!/bin/bash

set -e

git clone https://github.com/miracl/amcl.git
cd amcl/version3/c
git reset --hard e39201e3d34f4406530c103bb01f50fd84253b48

