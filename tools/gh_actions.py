#!/usr/bin/env python3
"""CLI minimo para GitHub Actions usando el token guardado por git.

Uso:
  python tools/gh_actions.py status [run_id] [workflow.yml]
  python tools/gh_actions.py watch <run_id> [segundos]
  python tools/gh_actions.py artifacts <run_id>
  python tools/gh_actions.py artifact <run_id> <nombre>
  python tools/gh_actions.py dispatch <workflow.yml>
"""
import argparse
import io
import json
import pathlib
import subprocess
import sys
import time
import urllib.request
import zipfile

REPO = "emanuelbertey/spacebit"
API = f"https://api.github.com/repos/{REPO}"


def token():
    r = subprocess.run(["git", "credential", "fill"],
                       input="protocol=https\nhost=github.com\n\n",
                       text=True, capture_output=True, check=True)
    for line in r.stdout.splitlines():
        if line.startswith("password="):
            return line.split("=", 1)[1]
    raise RuntimeError("no hay token en git credential")


TOKEN = token()


def api(path, method="GET", data=None):
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(
        f"{API}{path}" if path.startswith("/") else f"{API}/{path}",
        data=body, method=method,
        headers={"Authorization": f"Bearer {TOKEN}",
                 "Accept": "application/vnd.github+json",
                 "User-Agent": "gh-actions-cli"})
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read()
        return json.loads(raw) if raw else None


def jobs(run_id):
    return api(f"/actions/runs/{run_id}/jobs")["jobs"]


def show_run(r):
    print(f"run={r['id']} #{r.get('run_number')} {r['status']} / "
          f"{r.get('conclusion')} sha={r['head_sha'][:7]}")
    for j in jobs(r["id"]):
        print(f"  {j['name']}: {j['status']} / {j.get('conclusion')}")
        for s in j["steps"]:
            print(f"     {s['name']}: {s.get('conclusion') or '-'}")


def cmd_status(a):
    if a.run_id:
        show_run(api(f"/actions/runs/{a.run_id}"))
    else:
        path = f"/actions/workflows/{a.workflow}/runs" if a.workflow else "/actions/runs"
        for r in api(path + "?per_page=10")["workflow_runs"]:
            print(f"{r['id']} #{r['run_number']} {r['status']} / "
                  f"{r.get('conclusion')} {r['head_sha'][:7]}")


def cmd_watch(a):
    ultimo = None
    while True:
        r = api(f"/actions/runs/{a.run_id}")
        estado = (r["status"], r.get("conclusion"))
        if estado != ultimo:
            print(f"[{time.strftime('%H:%M:%S')}] {estado[0]} / {estado[1]}")
            ultimo = estado
        if r["status"] == "completed":
            for j in jobs(r["id"]):
                print(f"  {j['name']}: {j.get('conclusion')}")
            break
        time.sleep(a.interval)


def cmd_artifacts(a):
    for x in api(f"/actions/runs/{a.run_id}/artifacts")["artifacts"]:
        print(f"{x['name']} id={x['id']} size={x['size_in_bytes']}")


class AuthRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None:
            new.add_header("Authorization", f"Bearer {TOKEN}")
        return new


def cmd_artifact(a):
    arts = api(f"/actions/runs/{a.run_id}/artifacts")["artifacts"]
    hit = [x for x in arts if x["name"] == a.name]
    if not hit:
        print("no existe ese artefacto", file=sys.stderr)
        return 1
    req = urllib.request.Request(hit[0]["archive_download_url"],
                                 headers={"Authorization": f"Bearer {TOKEN}",
                                          "User-Agent": "gh-actions-cli"})
    import tempfile
    tmp = pathlib.Path(tempfile.gettempdir()) / f"artifact_{hit[0]['id']}.zip"
    subprocess.run(["curl.exe", "-sS", "-L", "-A", "gh-actions-cli",
                    "-H", f"Authorization: Bearer {TOKEN}",
                    "-o", str(tmp), hit[0]["archive_download_url"]], check=True)
    z = zipfile.ZipFile(tmp)
    for n in z.namelist():
        if a.file is None or n == a.file:
            print(f"=== {n} ===")
            print(z.read(n).decode("utf-8", "replace"))
    return 0


def cmd_dispatch(a):
    api(f"/actions/workflows/{a.workflow}/dispatches", "POST", {"ref": "main"})
    print("dispatch 204")


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")
    s = sub.add_parser("status")
    s.add_argument("run_id", nargs="?")
    s.add_argument("workflow", nargs="?", default=None)
    s.set_defaults(fn=cmd_status)
    w = sub.add_parser("watch")
    w.add_argument("run_id")
    w.add_argument("interval", nargs="?", type=int, default=30)
    w.set_defaults(fn=cmd_watch)
    ar = sub.add_parser("artifacts")
    ar.add_argument("run_id")
    ar.set_defaults(fn=cmd_artifacts)
    af = sub.add_parser("artifact")
    af.add_argument("run_id")
    af.add_argument("name")
    af.add_argument("--file", default=None)
    af.set_defaults(fn=cmd_artifact)
    d = sub.add_parser("dispatch")
    d.add_argument("workflow")
    d.set_defaults(fn=cmd_dispatch)
    a = p.parse_args()
    return a.fn(a) or 0


if __name__ == "__main__":
    sys.exit(main())
