# PFLOTRAN

PFLOTRAN is an open source, state-of-the-art massively parallel subsurface flow and reactive transport code. The code is developed under a GNU LGPL license allowing for third parties to interface proprietary software with the code; however, any modifications to the code itself must be documented and remain open source. PFLOTRAN is written in object-oriented, free-formatted Fortran 2003. The choice of Fortran over C/C++ was based primarily on the need to enlist and preserve tight collaboration with *experienced* domain scientists, without which PFLOTRAN's sophisticated process models would not exist.

## Official MCP server

The official Model Context Protocol (MCP) server for PFLOTRAN is maintained in a **separate** repository:

- **Canonical (Bitbucket):** https://bitbucket.org/pflotran/pflotran-mcp

`pflotran-mcp` is a Python control layer for agentic AI workflows, digital-twin training loops, verification, and output analysis. It wraps the **compiled** PFLOTRAN Fortran/MPI executable (via `PFLOTRAN_BIN` or a site module). PFLOTRAN itself remains a compiled application; the MCP package is not a Python reimplementation of the simulator.

Optional mirrors of the MCP project (for discovery only; development and issues should use Bitbucket) may appear on GitHub or GitLab under the `pflotran` organization. Those mirrors point back to the Bitbucket repository as the source of truth.

More detail: [README_MCP.md](README_MCP.md).

### Quick start (MCP)

```bash
git clone https://bitbucket.org/pflotran/pflotran-mcp.git
cd pflotran-mcp
python3 -m pip install -e ".[dev]"
export PFLOTRAN_BIN=/path/to/pflotran   # this tree's built executable
export PFLOTRAN_MCP_WORKSPACE=$HOME/pflotran-agent-runs
pflotran-mcp-server
```

## Related links

- PFLOTRAN website: https://www.pflotran.org/
- Documentation: https://documentation.pflotran.org/
- This repository (simulator): https://bitbucket.org/pflotran/pflotran
