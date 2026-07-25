# Official PFLOTRAN MCP Server

The Model Context Protocol (MCP) integration for PFLOTRAN lives in a **separate** repository.

## Canonical repository (source of truth)

**https://bitbucket.org/pflotran/pflotran-mcp**

Clone, install, open issues/PRs, and track releases there.

## What it is

`pflotran-mcp` is a Python **agent-facing control layer** around the compiled PFLOTRAN Fortran/MPI executable. It is **not** a Python port of PFLOTRAN.

```text
PFLOTRAN executable  →  pflotran-mcp  →  AI agent / workflow
```

Typical use: create cases, write input decks, run local or MPI jobs, inspect Tecplot/HDF5 output, and record provenance for digital-twin and verification workflows.

## This repository (PFLOTRAN)

This tree is the **simulator** (Fortran/MPI):

- Canonical: https://bitbucket.org/pflotran/pflotran

You still need a built `pflotran` binary (`PFLOTRAN_BIN`) for the MCP server to invoke.

## Mirrors

If you found `pflotran-mcp` on GitHub or GitLab, those are optional mirrors for discovery. Prefer Bitbucket for development.

| Host | MCP project |
|------|-------------|
| **Bitbucket (canonical)** | https://bitbucket.org/pflotran/pflotran-mcp |
| GitHub (mirror, if published) | https://github.com/pflotran/pflotran-mcp |
| GitLab (mirror, if published) | https://gitlab.com/pflotran/pflotran-mcp |

## Quick start

```bash
git clone https://bitbucket.org/pflotran/pflotran-mcp.git
cd pflotran-mcp
python3 -m pip install -e ".[dev]"
export PFLOTRAN_BIN=/path/to/pflotran
export PFLOTRAN_MCP_WORKSPACE=$HOME/pflotran-agent-runs
pflotran-mcp-server
```

See the MCP repository `README.md` and `docs/` for configuration, tools, and roadmap.
