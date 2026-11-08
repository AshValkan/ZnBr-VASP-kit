#!/bin/bash
# One-shot helper: build inputs for T1 (opt + opt_sol), preview the queue.
set -e
cd "$(dirname "$0")"
python3 bin/make_inputs.py --tasks opt,opt_sol
python3 bin/submit.py --tasks opt,opt_sol --dry-run
echo "--- check the list above, then run:"
echo "python3 bin/submit.py --tasks opt,opt_sol"
echo "python3 bin/status.py"
