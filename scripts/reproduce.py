"""Command-line entry point. No numerical replay runs by default or on import."""
import argparse
import ast
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
from common import ROOT, environment, frozen_check, output_path, require, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ["check", "export", "figures", "recompute"]:
        p = commands.add_parser(name)
        p.add_argument("--out", required=True, help="New output directory; existing paths are protected")
        if name == "export":
            p.add_argument("--records", help="Optional complete numerical replay directory")
        if name == "figures":
            p.add_argument("--analysis", required=True, help="Directory produced by export")
        if name == "recompute":
            p.add_argument("--confirm-frozen-replay", action="store_true")
            p.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.command == "recompute":
        from replay import recompute
        result = recompute(args.out, args.confirm_frozen_replay, args.resume)
    else:
        if args.command == "figures":
            frozen_check()
            analysis = output_path(args.analysis, create=False)
            require((analysis / "reporting_check.json").is_file(), "Run export first")
            out = output_path(args.out)
            result = subprocess.run(["Rscript", "--vanilla", str(ROOT / "scripts/plot_figures.R"),
                                     str(analysis), str(out / "figures")], check=True, capture_output=True, text=True)
            result = dict(accepted=True, figures=["Fig1.pdf", "Fig2.pdf"], new_gradient_updates=0)
            write(out / "figure_check.json", result)
        elif args.command == "export":
            from reporting import export
            out = output_path(args.out)
            result = export(out, output_path(args.records, create=False) if args.records else None)
        else:
            from records import validate
            result = validate()
            for path in [*ROOT.glob("scripts/*.py"), *ROOT.glob("code/**/*.py")]:
                ast.parse(path.read_text(), filename=str(path.relative_to(ROOT)))
            out = output_path(args.out)
            write(out / "saved_records_check.json", result)
        write(out / "environment.json", environment())
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
