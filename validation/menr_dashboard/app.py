import os
import streamlit as st
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import plotly.graph_objects as go

st.set_page_config(page_title="MENR Digital Twin", layout="wide")

BG = "#FFFFFF"
SURFACE = "#F4F4F2"
ACCENT = "#3FA796"
ACCENT_DIM = "#A8D9CE"
RED = "#C4543F"
TEXT = "#1A1D1F"
TEXT_DIM = "#5A5F63"
BORDER = "#E2E2DF"

st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@500;700&family=IBM+Plex+Mono:wght@400;500&family=Inter:wght@400;500&display=swap');

.stApp {{ background-color: {BG}; }}
html, body, [class*="css"] {{ font-family: 'Inter', sans-serif; }}

p, span, li, label, .stMarkdown, [data-testid="stCaptionContainer"] {{ color: {TEXT_DIM} !important; }}
h1, h2, h3, h4 {{ font-family: 'Space Grotesk', sans-serif; color: {TEXT} !important; }}
h2, h3 {{ color: {ACCENT} !important; }}

[data-testid="stMetric"] {{
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 16px 20px;
}}
[data-testid="stMetricValue"] {{ font-family: 'IBM Plex Mono', monospace; color: {ACCENT} !important; }}
[data-testid="stMetricLabel"] {{ color: {TEXT_DIM} !important; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }}

[data-testid="stTabs"] button p {{ color: {TEXT_DIM} !important; font-family: 'Space Grotesk', sans-serif; }}
[data-testid="stTabs"] button[aria-selected="true"] p {{ color: {ACCENT} !important; }}

[data-testid="stDataFrame"] {{ background-color: {SURFACE} !important; }}
[data-testid="stDataFrame"] * {{ color: {TEXT} !important; }}

[data-baseweb="select"] * {{ color: {TEXT} !important; }}
div[data-baseweb="tag"] {{ background-color: {ACCENT_DIM} !important; color: {TEXT} !important; }}

hr {{ border-color: {BORDER}; }}
</style>
""", unsafe_allow_html=True)

st.markdown(f"<p style='color:{ACCENT}; font-family:\"IBM Plex Mono\",monospace; font-size:13px; letter-spacing:0.08em; margin-bottom:4px;'>COMPUTATIONAL VALIDATION FRAMEWORK</p>", unsafe_allow_html=True)
st.markdown(f"<h1 style='margin-top:0; margin-bottom:4px; font-size:42px; color:{TEXT};'>MENR</h1>", unsafe_allow_html=True)
st.markdown(f"<p style='color:{TEXT_DIM}; font-size:16px; margin-top:0;'>Mechano-Epigenetic Nano-Rewriter — stiffness-triggered, two-wave epigenetic reprogramming, validated against 16,911 real TCGA patient-gene predictions</p>", unsafe_allow_html=True)
st.markdown(f"<hr style='border-color:{BORDER}; margin-top:20px; margin-bottom:28px;'>", unsafe_allow_html=True)

st.markdown("<h1 style='margin-bottom:0;'>MENR</h1>", unsafe_allow_html=True)
st.markdown("<p style='color:#5a6b5a; font-size:18px; margin-top:0;'>Mechano-Epigenetic Nano-Rewriter - Computational Validation Dashboard</p>", unsafe_allow_html=True)
st.markdown("---")

tab1, tab2, tab3, tab4 = st.tabs(["Model A - Stiffness Trigger", "Model B - Two-Wave Kinetics", "Digital Twin Lookup", "Sensitivity Sandbox"])

with tab1:
    st.subheader("Stiffness-triggered activation probability")

    cancers = {
        "Pancreatic (PAAD)":        {"E0": 3.5, "k": 0.8,  "color": "#2F6B4F"},
        "Triple-neg Breast (TNBC)": {"E0": 9.0, "k": 0.35, "color": "#5B8C5A"},
        "Hepatocellular (LIHC)":    {"E0": 5.5, "k": 0.6,  "color": "#8FAE6E"},
        "Lung Adeno (LUAD)":        {"E0": 5.5, "k": 0.6,  "color": "#B05C4A"},  # reused from LIHC - no solid lung MRE pair found (stated limitation)
    }

    col1, col2 = st.columns([1, 2])
    with col1:
        cancer_choice = st.selectbox("Cancer type", list(cancers.keys()))
        stiffness = st.slider("Tissue stiffness (kPa)", 0.0, 60.0, 20.0, 0.5)

        params = cancers[cancer_choice]
        prob = 1 / (1 + np.exp(-params["k"] * (stiffness - params["E0"])))
        st.metric("Activation probability", f"{prob*100:.1f}%")

        if stiffness <= 2.0:
            st.info("Soft / healthy tissue range")
        elif stiffness >= 10.0:
            st.success("Stiff / tumour tissue range")
        else:
            st.warning("Intermediate zone")

    with col2:
        E_values = np.linspace(0, 60, 300)
        fig = go.Figure()
        for name, p in cancers.items():
            P = 1 / (1 + np.exp(-p["k"] * (E_values - p["E0"])))
            fig.add_trace(go.Scatter(x=E_values, y=P, mode="lines", name=name,
                                      line=dict(width=3, color=p["color"])))
        fig.add_vline(x=stiffness, line_dash="dash", line_color="gray")
        fig.add_hline(y=0.5, line_dash="dot", line_color="lightgray")
        fig.add_vrect(x0=0.1, x1=2.0, fillcolor="#5B8C5A", opacity=0.12, line_width=0)
        fig.add_vrect(x0=10, x1=50, fillcolor="#B05C4A", opacity=0.08, line_width=0)
        fig.update_layout(
            xaxis_title="Tissue stiffness (kPa)", yaxis_title="Activation probability",
            template="plotly_white", paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF", height=450, legend=dict(orientation="h", y=-0.2)
        )
        st.plotly_chart(fig, use_container_width=True)

with tab2:
    st.subheader("Two-wave epigenetic reprogramming kinetics")

    col1, col2 = st.columns([1, 2])
    with col1:
        mechanism = st.radio("Silencing mechanism", ["Epigenetic (methylation)", "Genetic (deletion)"])
        m0 = st.slider("Starting methylation fraction", 0.0, 1.0, 0.65, 0.01)
        kme = st.slider("kme (DNMT3A rate)", 0.05, 0.45, 0.15, 0.01)
        kde = st.slider("kde (TET1 rate)", 0.05, 0.45, 0.40, 0.01)

    def DNMT3A_active(t, natural_level=0.12, inhibitor_strength=0.85, inhibitor_decay=0.4):
        inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
        return natural_level * (1 - inhibition)

    def TET1_active(t, delay=12, rise=0.3, clearance=0.08):
        if t < delay:
            return 0.0
        return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

    def methylation_ode(t, m, kme=kme, kde=kde):
        dnmt3a = DNMT3A_active(t)
        tet1 = TET1_active(t)
        return [kme * (1 - m[0]) * dnmt3a - kde * m[0] * tet1]

    if mechanism.startswith("Genetic"):
        t_eval = np.linspace(0, 48, 200)
        m_curve = np.full_like(t_eval, m0)
        final_m, pct_drop = m0, 0.0
        responsive = False
    else:
        sol = solve_ivp(methylation_ode, (0, 48), [m0], t_eval=np.linspace(0, 48, 200))
        t_eval, m_curve = sol.t, sol.y[0]
        final_m = m_curve[-1]
        pct_drop = (m0 - final_m) / m0 * 100 if m0 > 0 else 0
        responsive = pct_drop >= 40

    with col1:
        st.metric("Final methylation", f"{final_m:.3f}")
        st.metric("Percent drop", f"{pct_drop:.1f}%")
        st.metric("Predicted responsive (40%+ drop)", "Yes" if responsive else "No")

    with col2:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=t_eval, y=m_curve, mode="lines", line=dict(width=3, color=ACCENT),
                                  name="Methylation fraction"))
        fig.add_hline(y=m0*0.6, line_dash="dot", line_color="gray",
                      annotation_text="40% relative-drop threshold")
        fig.add_vline(x=12, line_dash="dash", line_color="lightgray", annotation_text="Wave 2 onset")
        fig.update_layout(xaxis_title="Time (hours)", yaxis_title="Promoter methylation fraction",
                           yaxis_range=[0, 1], template="plotly_white", paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF", height=450)
        st.plotly_chart(fig, use_container_width=True)

with tab3:
    st.subheader("Real patient-gene predictions (TCGA digital twin)")

    uploaded = st.file_uploader("Upload digital_twin_results_v2.csv", type="csv")
    df = None
    if uploaded is not None:
        df = pd.read_csv(uploaded)
    else:
        try:
            df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "digital-twin", "digital_twin_results_v2.csv"))
            st.caption("Loaded digital_twin_results_v2.csv from project folder.")
        except FileNotFoundError:
            st.warning("No results file found yet - upload digital_twin_results_v2.csv to explore real predictions.")

    if df is not None:
        col1, col2, col3 = st.columns(3)
        with col1:
            cancer_filter = st.multiselect("Cancer type", sorted(df["cancer_type"].unique()),
                                            default=sorted(df["cancer_type"].unique()))
        with col2:
            gene_filter = st.multiselect("Gene", sorted(df["gene"].unique()),
                                          default=sorted(df["gene"].unique()))
        with col3:
            mech_filter = st.multiselect("Mechanism", sorted(df["mechanism"].unique()),
                                          default=sorted(df["mechanism"].unique()))

        filtered = df[df["cancer_type"].isin(cancer_filter) &
                      df["gene"].isin(gene_filter) &
                      df["mechanism"].isin(mech_filter)]

        c1, c2, c3 = st.columns(3)
        c1.metric("Patients in filter", len(filtered))
        c2.metric("Responsive rate", f"{filtered['predicted_responsive'].mean()*100:.1f}%")
        c3.metric("Median starting methylation", f"{filtered['starting_methylation'].median():.3f}")

        fig = go.Figure()
        rates = filtered.groupby("cancer_type")["predicted_responsive"].mean() * 100
        fig.add_trace(go.Bar(x=rates.index, y=rates.values, marker_color=ACCENT))
        fig.update_layout(yaxis_title="Responsive rate (%)", template="plotly_white", paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF", height=400)
        st.plotly_chart(fig, use_container_width=True)

        st.dataframe(filtered.head(200), use_container_width=True)

with tab4:
    st.subheader("How fragile is the result")
    st.caption("Drag kme and kde to recreate the sensitivity analysis live.")

    col1, col2 = st.columns([1, 2])
    with col1:
        kme_mult = st.slider("kme multiplier", 0.5, 3.0, 1.0, 0.1, key="kme_s")
        kde_mult = st.slider("kde multiplier", 0.1, 2.0, 1.0, 0.1, key="kde_s")

    def DNMT3A_active2(t, natural_level=0.12, inhibitor_strength=0.85, inhibitor_decay=0.4):
        inhibition = inhibitor_strength * np.exp(-inhibitor_decay * t)
        return natural_level * (1 - inhibition)

    def TET1_active2(t, delay=12, rise=0.3, clearance=0.08):
        if t < delay:
            return 0.0
        return (1 - np.exp(-rise * (t - delay))) * np.exp(-clearance * (t - delay))

    def run_model(kme, kde, starts):
        results = []
        for m0_ in starts:
            def ode(t, m, kme=kme, kde=kde):
                dnmt = DNMT3A_active2(t)
                tet1 = TET1_active2(t)
                return [kme * (1 - m[0]) * dnmt - kde * m[0] * tet1]
            sol = solve_ivp(ode, (0, 48), [m0_], t_eval=[48])
            final = sol.y[0][-1]
            pct_drop = (m0_ - final) / m0_ * 100 if m0_ > 0 else 0
            results.append(pct_drop >= 40)
        return np.mean(results) * 100

    borderline = np.linspace(0.3, 0.65, 50)
    baseline_rate = run_model(0.15, 0.40, borderline)
    new_rate = run_model(0.15 * kme_mult, 0.40 * kde_mult, borderline)
    shift = new_rate - baseline_rate

    with col1:
        st.metric("Baseline responsive rate", f"{baseline_rate:.1f}%")
        st.metric("Adjusted responsive rate", f"{new_rate:.1f}%")
        st.metric("Shift from baseline", f"{shift:+.1f} pts")

    with col2:
        fig = go.Figure()
        fig.add_trace(go.Bar(x=["Baseline", "Adjusted"], y=[baseline_rate, new_rate],
                              marker_color=[ACCENT, RED if shift < 0 else ACCENT]))
        fig.update_layout(yaxis_title="Responsive rate among borderline patients (%)",
                           template="plotly_white", paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF", height=450)
        st.plotly_chart(fig, use_container_width=True)
