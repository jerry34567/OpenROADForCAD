#!/bin/bash

# apply STA patches
echo "Applying STA patches..."
cd src/sta
if git apply --check ../../my-sta-changes.patch 2>/dev/null; then
    git apply ../../my-sta-changes.patch
    echo "STA patches applied successfully"
elif git apply --reverse --check ../../my-sta-changes.patch 2>/dev/null; then
    echo "Patches already applied, skipping..."
else
    echo "Warning: Cannot apply patches cleanly, continuing anyway..."
fi
cd ../..

# compile
cd build 
make -j$(nproc)
make install
cd ..