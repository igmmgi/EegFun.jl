#!/bin/bash
# Post-build fix for VitePress theme imports
# This script must run AFTER julia docs/build.jl completes

THEME_FILE="docs/build/.vitepress/theme/index.ts"

if [ ! -f "$THEME_FILE" ]; then
    echo "Error: Theme file not found at $THEME_FILE"
    exit 1
fi

# Replace @/ imports with relative ./  imports
sed -i 's|from "@/VersionPicker.vue"|from "./VersionPicker.vue"|g' "$THEME_FILE"
sed -i "s|from '@/AuthorBadge.vue'|from './AuthorBadge.vue'|g" "$THEME_FILE"
sed -i "s|from '@/Authors.vue'|from './Authors.vue'|g" "$THEME_FILE"

# Copy Vue components from source to build theme directory
cp docs/src/.vitepress/theme/VersionPicker.vue docs/build/.vitepress/theme/
cp docs/src/.vitepress/theme/AuthorBadge.vue docs/build/.vitepress/theme/
cp docs/src/.vitepress/theme/Authors.vue docs/build/.vitepress/theme/

echo "✓ Fixed theme imports and copied Vue components"
