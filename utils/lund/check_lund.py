#!/usr/bin/env python3
import ROOT
import os
import sys
import glob
import math

# ----------------- Constants -----------------
PROTON_MASS = 0.938   # GeV
BEAM_ENERGY = 10.6    # GeV (massless e- along +z)

# ----------------- Globals -------------------
event_count = 0

# ----------------- Histograms ----------------
h2_q2_xbj = ROOT.TH2D("h2_q2_xbj",
                      "Q^{2} vs x_{Bj};x_{Bj};Q^{2} [GeV^{2}]",
                      100, 0.0, 1.5, 100, -1, 8.0)

h1_q2 = ROOT.TH1D("h1_q2",
                  "Q^{2} Distribution;Q^{2} [GeV^{2}];Counts",
                  100, -1, 4.0)

h1_q2_log = ROOT.TH1D("h1_q2_log",
                  "Q^{2} Distribution Log scale;Q^{2} [GeV^{2}];Counts",
                  100, -1, 8.0)


h1_xbj = ROOT.TH1D("h1_xbj",
                   "x_{Bj} Distribution;x_{Bj};Counts",
                   100, -1, 1.5)


# W and (Q^2 vs W)
h1_W = ROOT.TH1D("h1_W",
                 "W Distribution;W [GeV];Counts",
                 100, 0.5, 3.5)

h2_q2_W = ROOT.TH2D("h2_q2_W",
                    "Q^{2} vs W;W [GeV];Q^{2} [GeV^{2}]",
                    100, 0.5, 3.5, 100, -1, 8.0)

# --- New histograms ---
# Electron
h2_e_theta_p = ROOT.TH2D("h2_e_theta_p", "Electron #theta vs |p|;|p_{e}| [GeV];#theta_{e} [deg]", 100, 0, 7, 100, 0, 45)
h2_e_phi_p   = ROOT.TH2D("h2_e_phi_p", "Electron #phi vs |p|;|p_{e}| [GeV];#phi_{e} [deg]", 100, 0, 7, 100, -180, 180)
h2_e_theta_phi = ROOT.TH2D("h2_e_theta_phi", "Electron #theta vs #phi;#phi_{e} [deg];#theta_{e} [deg]", 100, -180, 180, 100, 0, 45)
h1_e_p = ROOT.TH1D("h1_e_p", "Electron |p| Distribution;|p_{e}| [GeV];Counts", 50, 0, 7)



# ---------------- Beam 4-vector ---------------
beam_vec = ROOT.TLorentzVector(0.0, 0.0, BEAM_ENERGY, BEAM_ENERGY)


def _parse_int(token, default=None):
    try:
        return int(token)
    except Exception:
        return default


def _parse_float(token, default=None):
    try:
        return float(token)
    except Exception:
        return default


def process_file(filepath):
    """
    Read a single LUND file and fill histograms.
    Assumes event header line followed by exactly N particle lines.
    Blank lines are ignored.
    """
    global event_count

    with open(filepath, "r") as f:
        # Strip whitespace and drop blanks so header & particles are contiguous
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    n_lines = len(lines)
    while i < n_lines:
        # --- Header ---
        header_tokens = lines[i].split()
        n_particles = _parse_int(header_tokens[0], default=None)
        if n_particles is None:
            i += 1
            continue

        # Incomplete event guard
        if i + 1 + n_particles > n_lines:
            break

        # ---- use header (DON'T SKIP IT) ----
        # Your header has 10 numbers: ... xB(5) W2(7) Q2(8) nu(9)
        if len(header_tokens) == 10:
            xB = _parse_float(header_tokens[5], default=-1.0)
            W2 = _parse_float(header_tokens[7], default=-1.0)
            Q2 = _parse_float(header_tokens[8], default=-1.0)
            nu = _parse_float(header_tokens[9], default=-1.0)
            
            #print("header")

            if W2:
                W = math.sqrt(W2)
                h1_W.Fill(W)
                h2_q2_W.Fill(W, Q2)
            if Q2:
                h2_q2_xbj.Fill(xB, Q2)
                h1_q2.Fill(Q2)
                h1_q2_log.Fill(Q2)
                h1_xbj.Fill(xB)

            event_count += 1  # count one processed event

        # --- Particle block (still only the N lines) ---
        event_lines = lines[i + 1: i + 1 + n_particles]
        i += 1 + n_particles  # advance cursor

        scattered_electron = None
        proton_momentum_mag = None

        for line in event_lines:
            parts = line.split()
            if len(parts) != 14:
                continue  # only particle lines here

            # LUND particle fields:
            # 0:index 1:lifetime 2:type 3:pid 4:parent 5:firstDau 6:px 7:py 8:pz 9:E 10:m 11:vx 12:vy 13:vz
            #print("particle")
            pid = _parse_int(parts[3], default=0)

            px = _parse_float(parts[6], default=0.0)
            py = _parse_float(parts[7], default=0.0)
            pz = _parse_float(parts[8], default=0.0)
            E  = _parse_float(parts[9], default=0.0)

            

            if pid == 11:
                p_mag = math.sqrt(px*px + py*py + pz*pz)
                theta = math.degrees(math.acos(pz / p_mag)) if p_mag > 0 else 0.0
                phi   = math.degrees(math.atan2(py, px))
                scattered_electron = ROOT.TLorentzVector(px, py, pz, E)
                h2_e_theta_p.Fill(p_mag, theta)
                h1_e_p.Fill(p_mag)          # make sure h1_e_p is defined!
                h2_e_phi_p.Fill(p_mag, phi)
                h2_e_theta_phi.Fill(phi, theta)


        # --- Kinematics from scattered electron ---
        #if scattered_electron:
        #    q = beam_vec - scattered_electron
        #    Q2 = -q.Mag2()
        #    nu = q.E()
        #    # --- Compute W from electron kinematics ---
        #    W2 = PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2
        #    if W2 > 0.0:
        #        W = math.sqrt(W2)
        #        h1_W.Fill(W)
        #        h2_q2_W.Fill(W, Q2)
        #if scattered_electron:
        #    W = math.sqrt(W2) if W2 > 0.0 else -1.0
        #    h1_W.Fill(W)
        #    h2_q2_W.Fill(W, Q2)
        #    
        #    if Q2 > 0.0 and nu > 0.0:
        #        xbj = Q2 / (2.0 * PROTON_MASS * nu)
        #        if xbj > 0.0:
        #            h2_q2_xbj.Fill(xbj, Q2)
        #            h1_q2.Fill(Q2)
        #            h1_xbj.Fill(xbj)
        #            if proton_momentum_mag is not None:
        #                h2_q2_pp.Fill(proton_momentum_mag, Q2)
        #            event_count += 1


def draw_with_text(hist, outname, text_x=0.2, text_y=0.85, logx=False, logy=False):
    c = ROOT.TCanvas("c_" + hist.GetName(), "", 900, 700)
    if logx: c.SetLogx()
    if logy:
        c.SetLogy()
        hist.SetMinimum(0.5)  # avoid zeros on log-y
    draw_opt = "COLZ" if hist.InheritsFrom("TH2") else ""
    hist.Draw(draw_opt)

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.04)
    lat.DrawLatex(text_x, text_y, f"N = {event_count}")

    c.SaveAs(outname)
    c.Close()



def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 check_lund.py /path/to/file.lund")
        print("  python3 check_lund.py '/path/to/dir/*.lund'")
        print("  python3 check_lund.py /path/to/dir")
        sys.exit(1)

    pattern = sys.argv[1]

    # Build a list of files: support directory, glob, or single file
    paths = []
    if os.path.isdir(pattern):
        paths = sorted(glob.glob(os.path.join(pattern, "*.lund")))
    else:
        # Pattern may be an explicit file or a glob
        matched = glob.glob(pattern)
        if matched:
            # Keep only files
            paths = sorted(p for p in matched if os.path.isfile(p))
        elif os.path.isfile(pattern):
            paths = [pattern]

    if not paths:
        print(f"No LUND files matched: {pattern}")
        sys.exit(1)

    print(f"Processing {len(paths)} file(s)")
    for p in paths:
        print(f"  -> {p}")
        process_file(p)
        
    preface = "incl_EG"

    # Draw/save plots once after all files are processed
    #draw_with_text(h2_q2_xbj, f"{preface}_Q2_vs_xbj.pdf")
    draw_with_text(h1_q2, f"{preface}_Q2_hist.pdf")
    #draw_with_text(h1_xbj, f"{preface}_xbj_hist.pdf")
    #draw_with_text(h1_W,     f"{preface}_W_hist.pdf")
    #draw_with_text(h2_q2_W,  f"{preface}_Q2_vs_W.pdf")
    draw_with_text(h1_e_p,  f"{preface}_electron momentum.pdf")
    #draw_with_text(h2_e_theta_p, f"{preface}_electron_theta_vs_p.pdf")
    #draw_with_text(h2_e_phi_p, f"{preface}_electron_phi_vs_p.pdf")
    draw_with_text(h1_q2_log, f"{preface}_Q2_hist_log.pdf", logy=True)
    
    #draw_with_text(h2_e_theta_phi, f"{preface}_electron_theta_vs_phi.pdf")
    ##draw_with_text(h2_p_theta_p, f"{preface}_proton_theta_vs_p.pdf")
    ##draw_with_text(h2_p_phi_p, f"{preface}_proton_phi_vs_p.pdf")
    ##draw_with_text(h2_p_theta_phi, f"{preface}_proton_theta_vs_phi.pdf")
    
    

    print("Done")


if __name__ == "__main__":
    # Ensure ROOT doesn’t pop interactive canvases
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    main()
