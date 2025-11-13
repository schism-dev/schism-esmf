---
title: Build Troubleshooting
summary: Known build issues and how to fix them
---

# Build Troubleshooting

For the complete list of known build issues with solutions and diagnostics, see the central troubleshooting document in the repository:

- View full troubleshooting: <https://github.com/schism-dev/schism-esmf/blob/HEAD/.github/BUILD_TROUBLESHOOTING.md>

Quick tips:

- Ensure ESMF is built with real MPI (not mpiuni)
- Verify SCHISM_BUILD_DIR points to a valid build with required libraries
- Use `cmake --build . -- VERBOSE=1` to diagnose link errors
- Run `.github/test-cmake-instructions.sh` to validate your environment
