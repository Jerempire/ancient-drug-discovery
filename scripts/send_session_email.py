"""Send ERAP2 session summary via Gmail SMTP (same pattern as ideas-hub)."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import os
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path

from dotenv import load_dotenv

# Load Gmail credentials from ideas-hub .env (shared)
load_dotenv(Path("C:/Users/jmj2z/Projects/intelligence/ideas-hub/.env"))
GMAIL_ADDRESS = os.getenv("GMAIL_ADDRESS", "")
GMAIL_APP_PASSWORD = os.getenv("GMAIL_APP_PASSWORD", "")
RECIPIENT = os.getenv("GMAIL_ADDRESS", "")  # send to self

SUBJECT = "ERAP2: The Full Story — From Black Death Genes to Peptide Drug Candidates"

HTML_BODY = """\
<!DOCTYPE html>
<html><body style="font-family: Georgia, serif; max-width: 720px; margin: 0 auto; line-height: 1.7; color: #222;">

<h1 style="border-bottom: 2px solid #333;">The ERAP2 Story: From Black Death Genes to Peptide Drug Candidates</h1>
<p style="color: #666; font-size: 14px;">A narrative of the ancient drug discovery project — March 2026</p>

<h2>The Thesis</h2>
<p>It started with a simple observation: the genes that helped our ancestors survive ancient plagues are the same genes that cancer uses to evade the immune system. ERAP2 — one of the strongest selection signals from the Black Death — trims peptides for immune presentation. If you could selectively block it, you might unlock immune responses against tumors that have learned to hide.</p>

<p>ERAP2 has a critical genetic variant: position 392 can be either lysine (K392, ancestral) or asparagine (N392, positively selected during the Black Death). K392 is associated with autoimmune risk. N392 processes peptides more efficiently. A drug that could distinguish between them would be textbook precision medicine.</p>

<h2>V1: The Naive Attempt</h2>
<p>We went straight at the active site — the obvious drug target. RFdiffusion designed 40 protein binders targeting the conserved catalytic residues (H370-393), ProteinMPNN gave us 320 sequences, and Boltz-2 said they all bound ERAP2 beautifully.</p>
<p>Then we ran the ERAP1 counter-screen. <strong>Every single binder preferred ERAP1 over ERAP2.</strong> Zero selectivity. The reason was painfully obvious in hindsight: 22 of 24 active site residues are identical between ERAP1 and ERAP2. We'd designed a perfect ERAP1 drug by accident.</p>
<p><em>Lesson: the counter-screen is the most valuable step in the pipeline.</em></p>

<h2>V2: Finding the Divergent Patch</h2>
<p>We pulled out BLOSUM62 alignments and found what V1 missed — two stretches in the substrate channel (residues 353-367 and 400-414) where ERAP2 actually differs from ERAP1. We redesigned 30 binders targeting these divergent hotspots.</p>
<p>This time, 3 candidates showed real selectivity. The lead — erap2v2_long_2 — hit ipTM 0.801 with a delta of +0.454 over ERAP1. Four constructs got locked for synthesis. Patent scaffold drafted with 12 claims.</p>

<h3>Then Physics Checked the ML</h3>
<p>PyRosetta interface analysis revealed the interfaces were physically real (dG -88 to -120 REU). But it also showed something Boltz-2 missed: <strong>IRAP cross-reactivity was real</strong>, not an ML artifact. 137 ERAP2 surface residues are divergent from ERAP1 but conserved with IRAP — a "trap" that protein binders fall into.</p>
<p>The alanine scan identified the fix: residues Y87 and Y89 contact IRAP but NOT ERAP2. A double mutant (Y87A_Y89A) killed IRAP binding while preserving ERAP2: ipTM 0.748 ERAP2, 0.186 IRAP. The cleanest panel we've ever produced.</p>
<p>Multi-sample reruns also killed Y4A — ANPEP cross-reactivity at 0.812 was hidden by single-sample noise. <strong>Single-sample Boltz-2 scores are not reliable for selectivity.</strong></p>

<h2>V3: The Channel Problem</h2>
<p>We tried to go further — variant-selective binders that could distinguish K392 from N392. But RFdiffusion fundamentally couldn't reach the substrate channel where position 392 sits. Three attempts, three failures. It designs convex surface-to-surface interfaces — it simply can't generate binders for concave channels.</p>

<h2>V4: Nature Already Solved This</h2>
<p>This is where the project pivoted to its most elegant approach.</p>
<p>Published biochemistry shows N392 excises hydrophobic P1 residues <strong>165 times faster</strong> than K392. The lysine's positive charge repels hydrophobic substrates. So if you flip it — put a negatively charged residue (glutamate or aspartate) at the P1 position of a peptide — it should form a salt bridge with K392's lysine and bind selectively.</p>
<p>Nature doesn't need a protein binder to tell these variants apart. It uses peptides. So that's what we did.</p>

<h3>DiffPepDock Results</h3>
<p>We designed a 28-peptide library with D-amino acid modifications (for protease resistance) and docked them via Boltz-2 + PyRosetta on Vast.ai. <strong>Salt bridge partially validated</strong>: 5/13 Glu/Asp peptides are K392-selective, but only 10+ mers — standard 9-mers universally prefer N392.</p>
<p>The biochemical explanation: longer peptides reach deeper into the K392 channel where the salt bridge actually forms. A 9-mer simply doesn't reach position 392.</p>

<h2>PepMLM: Independent Confirmation</h2>
<p>To test whether this was real, we brought in a completely independent tool: <strong>PepMLM</strong> (Nature Biotechnology 2025), a 650M-parameter model that generates peptide binders from target sequence alone.</p>
<p>PepMLM generated 48 scaffolds for the ERAP2 channel. It produced <strong>zero glutamate or aspartate at P1</strong> — the model strongly favors lysine (46%) and proline (19%). This means the salt bridge approach is genuinely novel: natural biology doesn't put negative charges at this position. Strong patent evidence.</p>
<p>We grafted our salt bridge chemistry (E/D at P1) onto PepMLM's ML-generated scaffolds, creating 32 hybrid peptides. <strong>7 out of 12 hybrids are K392-selective via Boltz-2</strong> — two independent scaffold sources, two different generation methods, same answer.</p>

<h3>Best Discoveries</h3>
<table style="border-collapse: collapse; width: 100%;">
<tr style="background: #f0f0f0;"><th style="padding: 8px; border: 1px solid #ccc;">Candidate</th><th style="padding: 8px; border: 1px solid #ccc;">Source</th><th style="padding: 8px; border: 1px solid #ccc;">K392 ipTM</th><th style="padding: 8px; border: 1px solid #ccc;">Selectivity</th></tr>
<tr><td style="padding: 8px; border: 1px solid #ccc;">pep_glu_long_01 (11-mer)</td><td style="padding: 8px; border: 1px solid #ccc;">Hand-designed</td><td style="padding: 8px; border: 1px solid #ccc;">0.774</td><td style="padding: 8px; border: 1px solid #ccc;">+0.097 (5-seed avg)</td></tr>
<tr><td style="padding: 8px; border: 1px solid #ccc;">hybrid_E_VKLLLL (10-mer)</td><td style="padding: 8px; border: 1px solid #ccc;">PepMLM hybrid</td><td style="padding: 8px; border: 1px solid #ccc;">0.798</td><td style="padding: 8px; border: 1px solid #ccc;">+0.075 (5-seed avg)</td></tr>
<tr><td style="padding: 8px; border: 1px solid #ccc;">hybrid_D_VKLLLL (10-mer)</td><td style="padding: 8px; border: 1px solid #ccc;">PepMLM hybrid</td><td style="padding: 8px; border: 1px solid #ccc;"><strong>0.812</strong></td><td style="padding: 8px; border: 1px solid #ccc;">+0.052</td></tr>
</table>

<h2>The Honest Check</h2>
<p>Multi-seed confirmation (5 seeds per candidate) showed the K392 selectivity <strong>direction is consistent</strong> but the <strong>magnitude lives in the noise band</strong>. More concerning: ALL peptides — including controls predicted to be N392-selective — lean K392 in Boltz-2. This could indicate a systematic scoring bias rather than real biology.</p>
<p>A bias check is running right now: 8 control peptides (random, scrambled, poly-Ala) + a P1 chemistry scan on the best scaffold. If random peptides also lean K392, the selectivity signal is weaker than it looks. If they don't, the signal is real.</p>

<h2>Where Things Stand</h2>
<h3>What's Real (High Confidence)</h3>
<ul>
<li>Glu at P1 is consistently the best salt bridge partner across all methods</li>
<li>Length >= 10 residues is required to reach position 392 in the channel</li>
<li>Two independent scaffold sources converge on the same selectivity direction</li>
<li>The salt bridge approach is genuinely novel (PepMLM never discovers it independently)</li>
<li>V2 binder Y87A_Y89A has the cleanest selectivity panel of any protein binder candidate</li>
</ul>

<h3>What's Uncertain</h3>
<ul>
<li>Whether K392 selectivity magnitude is real biology or Boltz-2 scoring bias</li>
<li>Whether N392-selective peptides exist (the evidence so far is shaky)</li>
<li>Whether D-amino acid modified peptides actually stay in the channel (kinetics unknown)</li>
</ul>

<h3>Two Viable Paths Forward</h3>
<ol>
<li><strong>V2 protein binder (Y87A_Y89A)</strong> — highest probability path. Ready for wet lab synthesis (~$3K). Blocks ERAP2 function for cancer immunotherapy.</li>
<li><strong>V4 peptides</strong> — higher upside, longer path. Variant-selective channel modulators for precision medicine. Pending bias check results before committing synthesis money (~$1,200).</li>
</ol>

<h2>Total Compute Cost: ~$7.70</h2>
<p>Across all Vast.ai runs from V1 through today's bias check. RTX 3090s and 4090s, $0.13-0.33/hr. The most expensive single run was V4 DiffPepDock at $0.60.</p>

<p style="color: #666; font-size: 12px; margin-top: 40px; border-top: 1px solid #ccc; padding-top: 10px;">Ancient Drug Discovery Pipeline — March 23, 2026<br>Generated by Claude Code</p>
</body></html>"""


def main():
    if not GMAIL_ADDRESS or not GMAIL_APP_PASSWORD:
        print("ERROR: Gmail credentials not found in ideas-hub .env")
        sys.exit(1)

    msg = MIMEMultipart("alternative")
    msg["Subject"] = SUBJECT
    msg["From"] = GMAIL_ADDRESS
    msg["To"] = RECIPIENT
    msg.attach(MIMEText(HTML_BODY, "html", "utf-8"))

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(GMAIL_ADDRESS, GMAIL_APP_PASSWORD)
        server.sendmail(GMAIL_ADDRESS, [RECIPIENT], msg.as_string())

    print(f"Email sent to {RECIPIENT}")


if __name__ == "__main__":
    main()
