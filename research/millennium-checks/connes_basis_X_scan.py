#!/usr/bin/env python3
"""Basis-size x X grid for connes_weil_form.Form (double precision, Legendre basis): |Mellin zero - true zero| for the first 5 zeros.
Separates the X effect from basis truncation (RH-6 confounded them). Writes connes_basis_X_scan.json.
"""
import json
from pathlib import Path

from connes_X_scan import scan

if __name__ == "__main__":
    res = {str(nb): scan((13, 17, 23, 29, 37, 47), nb) for nb in (24, 36, 48, 60, 72)}
    Path(__file__).with_name("connes_basis_X_scan.json").write_text(json.dumps(res, indent=1))
