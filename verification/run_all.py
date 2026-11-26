"""Run every verification script, store logs in logs/ and a summary in
logs/SUMMARY.txt.  Exit status is nonzero if any script fails.

The three programs in appendix/ are the certificate programs printed in the
paper's appendices (decoder, cone-identity verifier, Sturm-bound verifier);
the check_*.py scripts are independent checks written for this audit."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    ("appendix", "decode_certificate.py"),
    ("appendix", "verify_original_certificate.py"),
    ("appendix", "prove_square_scaling_29.py"),
    (".", "check_partitions.py"),
    (".", "check_reduction.py"),
    (".", "check_cone.py"),
    (".", "check_divisor.py"),
]


def main():
    logdir = os.path.join(HERE, "logs")
    os.makedirs(logdir, exist_ok=True)
    summary = []
    bad = 0
    for sub, s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, s], cwd=os.path.join(HERE, sub), capture_output=True, text=True)
        dt = time.time() - t0
        name = s[:-3]
        with open(os.path.join(logdir, name + ".log"), "w") as fh:
            fh.write(r.stdout)
            if r.stderr:
                fh.write("\n--- stderr ---\n" + r.stderr)
        last = [ln for ln in r.stdout.splitlines() if ln.startswith("--")]
        if last:
            line = last[-1]
        else:
            keys = ("Recovered", "Verified", "All ")
            hits = [ln for ln in r.stdout.splitlines() if ln.startswith(keys)]
            line = f"-- {name}: " + (hits[-1] if hits else "no summary")
        line += f"   [{dt:.1f}s]"
        if r.returncode != 0:
            bad += 1
            line += "   <-- FAILED"
        print(line, flush=True)
        summary.append(line)
    summary.append("ALL PASSED" if bad == 0 else f"{bad} SCRIPT(S) FAILED")
    with open(os.path.join(logdir, "SUMMARY.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print(summary[-1])
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
