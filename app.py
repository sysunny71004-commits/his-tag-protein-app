import re
import pandas as pd
import streamlit as st

st.set_page_config(page_title="Protein MW + Ni-NTA Buffer Designer", layout="wide")

# -------------------------
# 1) Protein MW calculator
# -------------------------
# Average residue masses (Da) for polypeptides (residue form, not free amino acids)
AA_RESIDUE_MW = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09,
    "C": 103.15, "E": 129.12, "Q": 128.13, "G": 57.05,
    "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
    "M": 131.19, "F": 147.18, "P": 97.12,  "S": 87.08,
    "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13
}

WATER = 18.015  # Da; lost per peptide bond when using free-AA masses. Here residue masses already assume that.
# NOTE: Because we use residue masses, we add terminal groups as H2O once (common convention).
# Many quick calculators do: sum(residue_masses) + H2O
# We'll do that for clarity.

def clean_sequence(seq: str) -> str:
    seq = seq.strip().upper()
    # Accept FASTA (remove header lines)
    seq = "\n".join([line for line in seq.splitlines() if not line.startswith(">")])
    seq = re.sub(r"[^A-Z]", "", seq)
    return seq

def calc_mw_da(seq: str) -> float:
    if not seq:
        return 0.0
    for aa in seq:
        if aa not in AA_RESIDUE_MW:
            raise ValueError(f"Invalid amino acid: {aa}")
    mw = sum(AA_RESIDUE_MW[aa] for aa in seq) + WATER
    return round(mw, 2)

def find_his_tag(seq: str):
    # Find longest run of H
    runs = [(m.start(), m.end()) for m in re.finditer(r"H{6,}", seq)]
    return runs  # list of (start, end)

# -------------------------
# 2) Buffer designer
# -------------------------
def ml_from_stock(final_vol_ml: float, final_conc, stock_conc) -> float:
    """C1V1 = C2V2. conc units must match. Returns mL."""
    return (final_conc * final_vol_ml) / stock_conc

def pct_from_stock(final_vol_ml: float, final_pct: float, stock_pct: float) -> float:
    """Percent v/v from a stock percent. Returns mL."""
    return (final_pct / stock_pct) * final_vol_ml

def make_buffer_table(components):
    df = pd.DataFrame(components, columns=["Component", "Stock", "Target (final)", "Volume (mL)"])
    df["Volume (mL)"] = df["Volume (mL)"].map(lambda x: round(x, 4))
    return df

def design_buffers(lysis_vol_ml, wash_vol_ml, elution_vol_ml,
                   use_tris_ph=8.0,
                   add_dtt=False, dtt_mM=1.0,
                   add_bme=False, bme_mM=5.0,
                   lysis_imid_mM=0.0):
    """
    Default design (QC-friendly typical):
    - Lysis: 50 mM Tris pH 8.0, 300 mM NaCl, 1% TritonX-100, PI 1x, Imidazole optional
    - Wash:  50 mM Tris pH 8.0, 300 mM NaCl, 20 mM Imidazole
    - Elution: 50 mM Tris pH 8.0, 300 mM NaCl, 250 mM Imidazole

    Available stocks:
    - Tris-HCl 1M pH 8.0 / 8.5
    - NaCl 5M
    - Triton X-100 10%
    - Protease inhibitor 100x
    - Imidazole 1M
    - DTT 1M
    - beta-me 14.3M
    """
    tris_stock = "1 M Tris-HCl (pH 8.0)" if use_tris_ph == 8.0 else "1 M Tris-HCl (pH 8.5)"
    tris_stock_conc_mM = 1000.0

    # ---- Lysis ----
    lysis = []
    v_tris = ml_from_stock(lysis_vol_ml, 50.0, tris_stock_conc_mM)  # 50 mM
    v_nacl = ml_from_stock(lysis_vol_ml, 300.0, 5000.0)            # 5M=5000 mM
    v_tx = pct_from_stock(lysis_vol_ml, 1.0, 10.0)                 # 10% stock to 1% final
    v_pi = lysis_vol_ml * (1.0/100.0)                              # 100x -> 1x
    v_imid = ml_from_stock(lysis_vol_ml, lysis_imid_mM, 1000.0) if lysis_imid_mM > 0 else 0.0

    lysis.append(["Tris-HCl (buffer base)", tris_stock, "50 mM", v_tris])
    lysis.append(["NaCl", "5 M NaCl", "300 mM", v_nacl])
    lysis.append(["Triton X-100", "10% Triton X-100", "1% (v/v)", v_tx])
    lysis.append(["Protease inhibitor", "100× PI", "1×", v_pi])
    if lysis_imid_mM > 0:
        lysis.append(["Imidazole (optional)", "1 M Imidazole", f"{lysis_imid_mM} mM", v_imid])

    v_dtt = 0.0
    if add_dtt and dtt_mM > 0:
        v_dtt = ml_from_stock(lysis_vol_ml, dtt_mM, 1000.0)
        lysis.append(["DTT (optional)", "1 M DTT", f"{dtt_mM} mM", v_dtt])

    v_bme = 0.0
    if add_bme and bme_mM > 0:
        # 14.3M stock => 14300 mM
        v_bme = ml_from_stock(lysis_vol_ml, bme_mM, 14300.0)
        lysis.append(["β-mercaptoethanol (optional)", "14.3 M β-ME", f"{bme_mM} mM", v_bme])

    used_lysis = v_tris + v_nacl + v_tx + v_pi + v_imid + v_dtt + v_bme
    lysis.append(["DW", "-", "to volume", max(0.0, lysis_vol_ml - used_lysis)])

    # ---- Wash ----
    wash = []
    v_tris_w = ml_from_stock(wash_vol_ml, 50.0, tris_stock_conc_mM)
    v_nacl_w = ml_from_stock(wash_vol_ml, 300.0, 5000.0)
    v_imid_w = ml_from_stock(wash_vol_ml, 20.0, 1000.0)

    wash.append(["Tris-HCl (buffer base)", tris_stock, "50 mM", v_tris_w])
    wash.append(["NaCl", "5 M NaCl", "300 mM", v_nacl_w])
    wash.append(["Imidazole", "1 M Imidazole", "20 mM", v_imid_w])

    used_wash = v_tris_w + v_nacl_w + v_imid_w
    wash.append(["DW", "-", "to volume", max(0.0, wash_vol_ml - used_wash)])

    # ---- Elution ----
    elution = []
    v_tris_e = ml_from_stock(elution_vol_ml, 50.0, tris_stock_conc_mM)
    v_nacl_e = ml_from_stock(elution_vol_ml, 300.0, 5000.0)
    v_imid_e = ml_from_stock(elution_vol_ml, 250.0, 1000.0)

    elution.append(["Tris-HCl (buffer base)", tris_stock, "50 mM", v_tris_e])
    elution.append(["NaCl", "5 M NaCl", "300 mM", v_nacl_e])
    elution.append(["Imidazole", "1 M Imidazole", "250 mM", v_imid_e])

    used_elution = v_tris_e + v_nacl_e + v_imid_e
    elution.append(["DW", "-", "to volume", max(0.0, elution_vol_ml - used_elution)])

    return make_buffer_table(lysis), make_buffer_table(wash), make_buffer_table(elution)

# -------------------------
# UI
# -------------------------
st.title("🧬 Protein MW Calculator + 🧪 Ni-NTA Buffer Designer")

col1, col2 = st.columns([1.1, 0.9], gap="large")

with col1:
    st.subheader("1) Protein sequence → Molecular weight")
    seq_in = st.text_area(
        "Amino-acid sequence (FASTA 가능)",
        height=220,
        placeholder="예) >protein\nMHHHHHHSSGLVPRGSHM...."
    )
    seq = clean_sequence(seq_in)

    cA, cB, cC = st.columns(3)
    with cA:
        show_his = st.checkbox("6×His-run 탐지 표시", value=True)
    with cB:
        assume_his_added = st.checkbox("N-말단 6×His가 추가된 단백질로 계산(+HHHHHH)", value=False)
    with cC:
        st.caption("※ residue MW 합 + H₂O(termini) 방식")

    if assume_his_added and seq:
        seq_calc = "HHHHHH" + seq
    else:
        seq_calc = seq

    if seq_calc:
        try:
            mw_da = calc_mw_da(seq_calc)
            st.metric("Molecular Weight", f"{mw_da:,.2f} Da", f"{mw_da/1000:,.3f} kDa")
            if show_his:
                runs = find_his_tag(seq_calc)
                if runs:
                    st.success(f"His-run (≥6) 탐지: {len(runs)}개  | 위치(0-index): {runs}")
                else:
                    st.info("His-run (≥6) 없음")
        except Exception as e:
            st.error(str(e))
    else:
        st.info("서열을 입력하면 분자량이 계산됩니다.")

with col2:
    st.subheader("2) Ni-NTA buffer designer (His-tag)")
    st.write("입력 조건(기본): **E. coli pellet 3 g**, **resin 100 µL (50% slurry)**")

    st.divider()
    st.markdown("**Buffer volume 설정**")
    lysis_vol = st.number_input("Lysis buffer final volume (mL)", min_value=5.0, max_value=500.0, value=30.0, step=5.0)
    wash_vol = st.number_input("Wash buffer final volume (mL)", min_value=2.0, max_value=500.0, value=10.0, step=2.0)
    elution_vol = st.number_input("Elution buffer final volume (mL)", min_value=2.0, max_value=500.0, value=10.0, step=2.0)

    st.divider()
    st.markdown("**옵션**")
    tris_ph = st.selectbox("Tris stock 선택", options=[8.0, 8.5], index=0, format_func=lambda x: f"1 M Tris-HCl pH {x}")
    lysis_imid = st.number_input("Lysis buffer imidazole (mM, optional)", min_value=0.0, max_value=50.0, value=0.0, step=5.0)

    add_dtt = st.checkbox("Add DTT", value=False)
    dtt_mM = st.number_input("DTT final (mM)", min_value=0.0, max_value=20.0, value=1.0, step=0.5, disabled=not add_dtt)

    add_bme = st.checkbox("Add β-ME", value=False)
    bme_mM = st.number_input("β-ME final (mM)", min_value=0.0, max_value=50.0, value=5.0, step=1.0, disabled=not add_bme)

    st.divider()
    run = st.button("✅ Calculate buffers", use_container_width=True)

if run:
    lysis_df, wash_df, elution_df = design_buffers(
        lysis_vol_ml=lysis_vol,
        wash_vol_ml=wash_vol,
        elution_vol_ml=elution_vol,
        use_tris_ph=tris_ph,
        add_dtt=add_dtt,
        dtt_mM=dtt_mM,
        add_bme=add_bme,
        bme_mM=bme_mM,
        lysis_imid_mM=lysis_imid
    )

    st.subheader("Results")
    tab1, tab2, tab3 = st.tabs(["Lysis buffer", "Wash buffer", "Elution buffer"])

    with tab1:
        st.dataframe(lysis_df, use_container_width=True)
    with tab2:
        st.dataframe(wash_df, use_container_width=True)
    with tab3:
        st.dataframe(elution_df, use_container_width=True)

    # Provide a code snippet that reproduces tables (useful for reports)
    st.subheader("Python snippet (report용 복붙)")
    snippet = f"""# Volumes (mL)
lysis_vol = {lysis_vol}
wash_vol = {wash_vol}
elution_vol = {elution_vol}

# Using Tris pH {tris_ph}
# Lysis imidazole = {lysis_imid} mM
# DTT: {'on' if add_dtt else 'off'} ({dtt_mM} mM)
# β-ME: {'on' if add_bme else 'off'} ({bme_mM} mM)
"""
    st.code(snippet, language="python")

st.caption("Made for His-tag Ni-NTA workflows (E. coli). You can extend to phosphate/HEPES bases if needed.")
