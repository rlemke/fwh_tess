"""Run a tess workflow in memory: the real FFL, the real evaluator, no MongoDB.

Every step waiting on an event facet is handed to ``handlers`` (qualified facet
name -> handler) and the workflow resumed, until it finishes. The runner does
the same dispatch for real; this is how the FFL is tested without one.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import facetwork
from facetwork import emit_dict, parse
from facetwork.runtime import Evaluator, ExecutionStatus, MemoryStore, Telemetry

FFL = Path(__file__).resolve().parents[1] / "src" / "tesstransit" / "ffl" / "tess.ffl"
HTTP = Path(facetwork.__file__).parent / "ffl" / "http.ffl"


def program() -> dict[str, Any]:
    return emit_dict(parse(HTTP.read_text() + "\n" + FFL.read_text()))


def _find(prog: dict[str, Any], qualified: str) -> dict[str, Any]:
    ns, name = qualified.rsplit(".", 1)
    for d in prog["declarations"]:
        if d.get("type") == "Namespace" and d["name"] == ns:
            for x in d["declarations"]:
                if x.get("type") == "WorkflowDecl" and x["name"] == name:
                    return x
    raise KeyError(qualified)


def run(
    workflow: str, inputs: dict[str, Any], handlers: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """``(workflow returns, facet names in the order they ran)``."""
    prog = program()
    wf = _find(prog, workflow)
    store = MemoryStore()
    ev = Evaluator(persistence=store, telemetry=Telemetry(enabled=False))
    res = ev.execute(wf, inputs=inputs, program_ast=prog)
    ran: list[str] = []
    while res.status == ExecutionStatus.PAUSED:
        blocked = [
            s
            for s in store._steps.values()
            if s.workflow_id == res.workflow_id and s.state == "state.EventTransmit"
        ]
        if not blocked:
            break
        for s in blocked:
            params = {k: v.value for k, v in s.attributes.params.items()}
            ran.append(s.facet_name)
            ev.continue_step(s.id, handlers[s.facet_name](params))
        res = ev.resume(res.workflow_id, wf, prog)
    if res.status != ExecutionStatus.COMPLETED or not res.success:
        raise RuntimeError(
            f"{workflow} ended {res.status}: {getattr(res, 'error', '')}"
        )
    return dict(res.outputs), ran
