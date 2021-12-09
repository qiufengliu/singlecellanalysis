#!/bin/bash

# add fetch and run batch helper script
chmod a+x /opt/ecs-additions/fetch_and_run.sh
cp /opt/ecs-additions/fetch_and_run.sh /usr/local/bin

# add awscli-shim
mv /opt/aws-cli/bin /opt/aws-cli/dist
chmod a+x /opt/ecs-additions/awscli-shim.sh
mkdir /opt/aws-cli/bin
cp /opt/ecs-additions/awscli-shim.sh /opt/aws-cli/bin/aws

# Remove current symlink
rm -f /usr/local/aws-cli/v2/current/bin/aws
cp /opt/ecs-additions/awscli-shim.sh /usr/local/aws-cli/v2/current/bin/aws # Used in Cromwell

# ensure that /usr/bin/aws points to the non-shimmed version
ln -sf /usr/local/aws-cli/v2/current/dist/aws /usr/bin/aws
