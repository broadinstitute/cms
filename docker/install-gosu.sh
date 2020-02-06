#!/bin/bash
set -exu -o pipefail

GOSU_VERSION=1.11

dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')"
wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch"
wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch.asc"

# verify the signature
export GNUPGHOME="$(mktemp -d)"
GPG_KEY=B42F6819007F00F88E364FD4036A9C25BF357DD4
(    gpg --keyserver ha.pool.sks-keyservers.net --recv-keys "$GPG_KEY" \
  || gpg --keyserver pool.sks-keyservers.net --recv-keys "$GPG_KEY" \
  || gpg --keyserver pgp.mit.edu --recv-keys "$GPG_KEY" \
  || gpg --keyserver keyserver.pgp.com --recv-keys "$GPG_KEY" )

gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu
rm -r "$GNUPGHOME" /usr/local/bin/gosu.asc
chmod +x /usr/local/bin/gosu

# verify that the binary works
gosu nobody true
