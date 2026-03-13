# Documentation Build Guide

This document explains how to build and work with the SCHISM-ESMF documentation.

## Overview

The project documentation is built using [MkDocs](https://www.mkdocs.org/) with the [Material theme](https://squidfunk.github.io/mkdocs-material/). The documentation sources are in Markdown format located in the `doc/` directory.

## Documentation Structure

```
schism-esmf/
├── doc/                          # Documentation source files
│   ├── index.md                  # Homepage
│   ├── installation.md           # Installation guide
│   ├── cmake-architecture.md     # CMake build system details
│   ├── running-examples.md       # How to run examples
│   ├── legal.md                  # License and legal information
│   ├── logo.png                  # Project logo
│   └── CMakeLists.txt            # CMake build rules for docs
├── mkdocs.yml                    # MkDocs configuration
├── Readme.md                     # Main project README (linked in docs)
├── QUICKSTART.md                 # Quick start guide (linked in docs)
└── .github/
    └── BUILD_TROUBLESHOOTING.md  # Build troubleshooting (linked in docs)
```

## Prerequisites

### Required

- **Python 3.8+**
- **MkDocs**: `pip install mkdocs`
- **MkDocs Material theme**: `pip install mkdocs-material`

### Installation Options

#### Using pip (recommended)

```bash
pip install mkdocs mkdocs-material
```

#### Using conda/mamba

```bash
conda install -c conda-forge mkdocs mkdocs-material
```

#### Using Poetry (for developers)

```bash
poetry install
```

## Building Documentation

### Using CMake (Integrated Build)

The documentation build is integrated into the CMake build system.

#### Enable Documentation Build

```bash
mkdir -p build && cd build
cmake .. -DBUILD_DOCS=ON
```

#### Build HTML Documentation

```bash
make docs
```

The built documentation will be in `build/docs/`.

#### Serve Documentation Locally

```bash
make docs-serve
```

This starts a local web server at http://127.0.0.1:8000 with live reloading.

#### Clean Documentation

```bash
make docs-clean
```

### Using MkDocs Directly

You can also build documentation directly with MkDocs without CMake:

#### Build

```bash
mkdocs build
```

Output goes to `site/` directory.

#### Serve Locally

```bash
mkdocs serve
```

Navigate to http://127.0.0.1:8000 to view.

#### Build for Deployment

```bash
mkdocs build --strict
```

The `--strict` flag treats warnings as errors.

## CMake Integration Details

### Build Option

The documentation build is controlled by the `BUILD_DOCS` CMake option:

```cmake
option(BUILD_DOCS "Build HTML documentation with MkDocs" OFF)
```

It's **OFF by default** to avoid requiring Python/MkDocs for users who only want to build the software.

### Docs-only Configure (no Fortran/MPI required)

If your Fortran/MPI toolchain isn't available, you can configure a documentation-only build:

```bash
mkdir -p build && cd build
cmake .. -DDOCS_ONLY=ON  # Implicitly sets -DBUILD_DOCS=ON

# Build user guide (MkDocs)
cmake --build . --target docs

# Build API reference (FORD - requires ford installed)
cmake --build . --target docs-api

# Build both
cmake --build . --target docs-all
```

For details on writing FORD-style documentation in the Fortran source code, see the [API Documentation Guide](api-documentation-guide.md).

### CMake Targets

| Target | Description |
|--------|-------------|
| `docs` | Build HTML user guide (MkDocs) |
| `docs-api` | Build API reference from Fortran source (FORD) |
| `docs-all` | Build both user guide and API reference |
| `docs-serve` | Start local documentation server |
| `docs-clean` | Remove built documentation |

### Installation

When you run `make install`, the documentation is installed to:

- **HTML documentation**: `${CMAKE_INSTALL_PREFIX}/share/doc/SCHISM_ESMF_Interface/html/`
- **Markdown sources**: `${CMAKE_INSTALL_PREFIX}/share/doc/SCHISM_ESMF_Interface/markdown/`
- **Top-level docs**: `${CMAKE_INSTALL_PREFIX}/share/doc/SCHISM_ESMF_Interface/` (Readme.md, QUICKSTART.md)

Example:

```bash
cmake .. -DBUILD_DOCS=ON -DCMAKE_INSTALL_PREFIX=/usr/local
make docs
sudo make install
```

## Continuous Integration

### GitHub Actions Example

```yaml
name: Build Documentation

on:
  push:
    branches: [main, develop]
  pull_request:

jobs:
  docs:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      
      - name: Install dependencies
        run: |
          pip install mkdocs mkdocs-material
      
      - name: Build documentation
        run: mkdocs build --strict
      
      - name: Deploy to GitHub Pages
        if: github.ref == 'refs/heads/main'
        run: mkdocs gh-deploy --force
```

### GitHub Pages Deployment

To deploy documentation to GitHub Pages:

```bash
mkdocs gh-deploy
```

This builds and pushes the documentation to the `gh-pages` branch.

## Writing Documentation

### File Format

All documentation files are Markdown (`.md`) with YAML front matter:

```markdown
---
title: Page Title
summary: Brief description
---

# Main Heading

Content here...
```

### Code Blocks

Use fenced code blocks with language specification:

````markdown
```python
def hello():
    print("Hello, SCHISM!")
```
````

### Admonitions

Use Material theme admonitions for callouts:

```markdown
!!! note "Important Note"
    This is an important note.

!!! warning
    This is a warning.

!!! tip
    This is a helpful tip.
```

### Internal Links

Link to other documentation pages:

```markdown
See [Installation Guide](installation.md) for details.
```

### External Links

Link to GitHub files:

```markdown
See [BUILD_TROUBLESHOOTING.md](../.github/BUILD_TROUBLESHOOTING.md).
```

## Troubleshooting

### Python Not Found

**Problem**: `Python3 not found` during CMake configure.

**Solution**: Install Python 3.8+ and ensure it's in your PATH.

```bash
# macOS
brew install python@3.11

# Ubuntu/Debian
sudo apt install python3 python3-pip

# Windows
# Download from python.org
```

### MkDocs Not Found

**Problem**: `MkDocs not found` warning during CMake configure.

**Solution**: Install MkDocs:

```bash
pip install mkdocs mkdocs-material
```

### Material Theme Not Found

**Problem**: `MkDocs Material theme not found` warning.

**Solution**: Install the Material theme:

```bash
pip install mkdocs-material
```

### Build Fails with Import Errors

**Problem**: Documentation build fails with Python import errors.

**Solution**: Check that all required packages are installed:

```bash
python3 -c "import mkdocs; print(mkdocs.__version__)"
python3 -c "import material; print(material.__version__)"
```

If either fails, reinstall:

```bash
pip install --upgrade mkdocs mkdocs-material
```

### Documentation Doesn't Update

**Problem**: Changes to markdown files don't appear in built docs.

**Solution**: 
1. Clean the build: `make docs-clean` or `rm -rf site/`
2. Rebuild: `make docs` or `mkdocs build`

When using `mkdocs serve`, changes should auto-reload. If not, restart the server.

## Configuration Reference

### mkdocs.yml Structure

```yaml
site_name: Project Name           # Site title
site_url: https://example.com     # Base URL
docs_dir: ./doc                   # Source directory
site_dir: ./site                  # Output directory
theme:
  name: material                  # Theme name
  features: [...]                 # Theme features
nav:                              # Navigation structure
  - Home: index.md
  - Guide: guide.md
markdown_extensions: [...]        # Markdown processing extensions
plugins: [...]                    # MkDocs plugins
```

### Material Theme Features

Enabled features in `mkdocs.yml`:

- `toc.integrate` - Integrate table of contents into navigation
- `toc.follow` - Follow active item in ToC
- `content.code.annotate` - Code annotation support
- `navigation.tabs` - Top-level navigation tabs
- `navigation.sections` - Group pages into sections
- `navigation.expand` - Expand sections by default
- `navigation.top` - "Back to top" button
- `search.suggest` - Search suggestions
- `search.highlight` - Highlight search terms

### Markdown Extensions

Key extensions enabled:

- `admonition` - Note/warning boxes
- `codehilite` - Syntax highlighting
- `pymdownx.superfences` - Nested code blocks
- `pymdownx.tabbed` - Tabbed content
- `pymdownx.details` - Collapsible sections

## Best Practices

1. **Keep files focused**: Each page should cover one topic
2. **Use clear headings**: Hierarchical structure with H1, H2, H3
3. **Test locally**: Always run `mkdocs serve` to preview changes
4. **Link generously**: Cross-link related content
5. **Use code examples**: Show, don't just tell
6. **Add screenshots**: Visual aids help comprehension (when applicable)
7. **Update navigation**: Add new pages to `mkdocs.yml` nav section
8. **Follow naming conventions**: Use lowercase-with-hyphens for filenames

## Maintenance

### Updating Dependencies

Check for MkDocs updates:

```bash
pip list --outdated | grep mkdocs
```

Update:

```bash
pip install --upgrade mkdocs mkdocs-material
```

### Checking for Broken Links

Install link checker:

```bash
pip install linkchecker
```

Check built docs:

```bash
mkdocs build
linkchecker site/
```

### Documentation Review Checklist

Before committing documentation changes:

- [ ] Spell-check content
- [ ] Verify all links work
- [ ] Build with `--strict` flag passes
- [ ] Preview with `mkdocs serve` looks correct
- [ ] Code examples are tested and correct
- [ ] Navigation updated if new pages added
- [ ] Front matter (title, summary) is accurate

## Resources

- [MkDocs Documentation](https://www.mkdocs.org/)
- [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
- [Markdown Guide](https://www.markdownguide.org/)
- [PyMdown Extensions](https://facelessuser.github.io/pymdown-extensions/)

## Getting Help

- **MkDocs Issues**: https://github.com/mkdocs/mkdocs/issues
- **Material Theme Issues**: https://github.com/squidfunk/mkdocs-material/issues
- **Project Docs Issues**: https://github.com/schism-dev/schism-esmf/issues

---

**Last Updated**: November 2025  
**Maintainer**: Carsten Lemmen <carsten.lemmen@hereon.de>
