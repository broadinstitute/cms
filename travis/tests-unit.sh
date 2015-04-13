#!/bin/bash
set -e

 nosetests -v --with-xunit --with-coverage --nocapture \
     --cover-inclusive --cover-branches --cover-tests \
     --cover-package selection \
     -w cms/test/unit/