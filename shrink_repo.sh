#!/bin/bash
set -e

echo "1. Committing current changes..."
git add .
git commit -m "Update docs and benchmarking" || echo "No changes to commit"

echo "2. Backing up binary files to a tarball..."
tar -cf binary_backup.tar resources/data

echo "3. Running git-filter-repo to purge binary extensions from all history..."
# git-filter-repo requires a clean working tree, which we now have.
./git-filter-repo --force \
  --invert-paths \
  --path-glob '*.bdf' \
  --path-glob '*.eeg' \
  --path-glob '*.vhdr' \
  --path-glob '*.vmrk' \
  --path-glob '*.edf' \
  --path-glob '*.fdt' \
  --path-glob '*.set' \
  --path-glob '*.mat' \
  --path-glob '*.fif' \
  --path-glob '*.jld2' \
  --path-glob '*.xdf'

echo "4. Restoring binary files from backup..."
tar -xf binary_backup.tar
rm binary_backup.tar

echo "5. Committing binary files as brand new..."
git add resources/data
git commit -m "Add final binary datasets"

echo "6. Cleaning up..."
rm git-filter-repo
rm shrink_repo.sh

echo "Done! Run 'du -sh .git' to see your new repository size!"
