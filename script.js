let sfdChart = null;
let bmdChart = null;
let lastData = null;

const fmt = (n, d = 2) =>
  Number.isFinite(n)
    ? n.toLocaleString("id-ID", { minimumFractionDigits: d, maximumFractionDigits: d })
    : "-";

// Faktor beta1 sesuai SNI 2847:2019 Tabel 22.2.2.4.3
const getBeta1 = (fc) =>
  fc <= 28 ? 0.85 : Math.max(0.65, 0.85 - (0.05 * (fc - 28)) / 7);

// Tinggi minimum balok sesuai SNI 2847:2019 Tabel 9.3.1.1
function getHminInfo(Ln_mm, support) {
  if (support === "fixed-fixed") return { hMin: Ln_mm / 21, formula: "Ln / 21 (Kedua ujung menerus)" };
  if (support === "fixed-pinned") return { hMin: Ln_mm / 18.5, formula: "Ln / 18,5 (Satu ujung menerus)" };
  if (support === "pinned-pinned") return { hMin: Ln_mm / 16, formula: "Ln / 16 (Tumpuan sederhana)" };
  return { hMin: Ln_mm / 8, formula: "Ln / 8 (Kantilever)" };
}

// Analisis gaya dalam berdasarkan kondisi tumpuan
function calcForces(Ln, qu, support) {
  let MuNeg = 0, MuPos = 0, VuMax = 0;
  if (support === "fixed-fixed") {
    MuNeg = (qu * Ln * Ln) / 12;
    MuPos = (qu * Ln * Ln) / 24;
    VuMax = (qu * Ln) / 2;
  } else if (support === "fixed-pinned") {
    MuNeg = (qu * Ln * Ln) / 8;
    MuPos = (9 * qu * Ln * Ln) / 128;
    VuMax = (5 * qu * Ln) / 8;
  } else if (support === "pinned-pinned") {
    MuNeg = 0;
    MuPos = (qu * Ln * Ln) / 8;
    VuMax = (qu * Ln) / 2;
  } else {
    MuNeg = (qu * Ln * Ln) / 2;
    MuPos = 0;
    VuMax = qu * Ln;
  }

  const x = [], V = [], M = [];
  const steps = 30;
  for (let i = 0; i <= steps; i++) {
    const xi = (Ln * i) / steps;
    x.push(xi.toFixed(1));
    V.push(Number((VuMax - qu * xi).toFixed(2)));
    M.push(Number((-MuNeg + VuMax * xi - 0.5 * qu * xi * xi).toFixed(2)));
  }
  return { MuNeg, MuPos, VuMax, x, V, M };
}

// Estimasi otomatis dimensi penampang (b x h) sesuai SNI 2847:2019
function estimateDimensions(Ln, support, qDL, qLL, fc, gammaC, fy, cc, ds, db) {
  const { hMin, formula } = getHminInfo(Ln * 1000, support);

  let h = Math.max(300, Math.ceil((hMin * 1.45) / 50) * 50);
  let b = Math.max(200, Math.round((0.55 * h) / 50) * 50);

  const beta1 = getBeta1(fc);
  const rhoTarget = Math.min(0.012, (3 / 8) * beta1 * ((0.85 * fc) / fy) * 0.65);

  for (let iter = 0; iter < 15; iter++) {
    const d = h - cc - ds - db / 2;
    const qSw = (b / 1000) * (h / 1000) * gammaC;
    const qu = 1.2 * (qDL + qSw) + 1.6 * qLL;
    const forces = calcForces(Ln, qu, support);
    const MuMax = Math.max(forces.MuNeg, forces.MuPos);
    const Rn = (MuMax * 1e6) / (0.9 * b * d * d);
    const m = fy / (0.85 * fc);
    const disc = 1 - (2 * m * Rn) / fy;
    const rhoReq = disc > 0 ? (1 / m) * (1 - Math.sqrt(disc)) : 1;

    const Vu_d = Math.max(0, forces.VuMax - qu * (d / 1000));
    const Vc = (0.17 * Math.sqrt(fc) * b * d) / 1000;
    const VsMax = (0.66 * Math.sqrt(fc) * b * d) / 1000;
    const VsReq = Vu_d > 0.75 * Vc ? (Vu_d - 0.75 * Vc) / 0.75 : 0;

    if (rhoReq <= rhoTarget && VsReq <= VsMax && h >= hMin) {
      return { b, h, d, hMin, formula };
    }

    h += 50;
    b = Math.max(200, Math.round((0.55 * h) / 50) * 50);
  }

  return { b, h, d: h - cc - ds - db / 2, hMin, formula };
}

function designFlexure(Mu_kNm, b, d, fc, fy, db) {
  const beta1 = getBeta1(fc);
  const m = fy / (0.85 * fc);
  const rhoMin = Math.max((0.25 * Math.sqrt(fc)) / fy, 1.4 / fy);
  const As1 = (Math.PI / 4) * db * db;

  let AsReq = 2 * As1;
  if (Mu_kNm > 0) {
    const Rn = (Mu_kNm * 1e6) / (0.9 * b * d * d);
    const disc = Math.max(0.05, 1 - (2 * m * Rn) / fy);
    const rho = (1 / m) * (1 - Math.sqrt(disc));
    AsReq = Math.max(rho, rhoMin) * b * d;
  }

  let n = Math.max(2, Math.ceil(AsReq / As1));
  let AsProv = n * As1;
  let a = (AsProv * fy) / (0.85 * fc * b);
  let c = a / beta1;
  let et = 0.003 * ((d - c) / c);
  let phi = et >= 0.005 ? 0.9 : Math.max(0.65, 0.65 + 0.25 * ((et - fy / 200000) / (0.005 - fy / 200000)));
  let phiMn = (phi * AsProv * fy * (d - a / 2)) / 1e6;

  while (phiMn < Mu_kNm && n < 12) {
    n++;
    AsProv = n * As1;
    a = (AsProv * fy) / (0.85 * fc * b);
    c = a / beta1;
    et = 0.003 * ((d - c) / c);
    phi = et >= 0.005 ? 0.9 : Math.max(0.65, 0.65 + 0.25 * ((et - fy / 200000) / (0.005 - fy / 200000)));
    phiMn = (phi * AsProv * fy * (d - a / 2)) / 1e6;
  }

  return { Mu_kNm, AsReq, n, AsProv, phiMn, et, phi, ok: phiMn >= Mu_kNm && et >= 0.004 };
}

/**
 * Desain Tulangan Geser per Zona (SNI 2847:2019 Bab 22, Tabel 9.7.6.2.2 & Detailing Zona Tumpuan)
 * zone: "tumpuan" atau "lapangan"
 */
function designShear(Vu_kN, b, d, fc, fyt, ds, zone) {
  const phi = 0.75;
  const Vc = (0.17 * Math.sqrt(fc) * b * d) / 1000;
  const phiVc = phi * Vc;
  const Av = 2 * (Math.PI / 4) * ds * ds;

  const VsLimit1 = (0.33 * Math.sqrt(fc) * b * d) / 1000;
  const VsReq = Vu_kN > phiVc ? (Vu_kN - phiVc) / phi : 0;

  // Batas jarak maksimum berdasarkan SNI 2847:2019 Tabel 9.7.6.2.2 & Pasal 9.6.3.3
  let sMaxCode = VsReq <= VsLimit1 ? Math.min(d / 2, 600) : Math.min(d / 4, 300);
  const sMaxMinReinf = (Av * fyt) / (Math.max(0.062 * Math.sqrt(fc), 0.35) * b);
  sMaxCode = Math.min(sMaxCode, sMaxMinReinf);

  // Pada zona tumpuan (daerah kritis geser & detailing ujung balok), diterapkan batas praktis maks 150 mm
  // Pada zona lapangan, mengikuti batas sMaxCode (misal d/2 = 195,25 mm -> 175 mm)
  const sMaxEffective = zone === "tumpuan" ? Math.min(sMaxCode, 150) : sMaxCode;

  let s = Math.max(50, Math.floor(sMaxEffective / 25) * 25);
  if (VsReq > 0) {
    const sReq = (Av * fyt * d) / (VsReq * 1000);
    s = Math.max(50, Math.min(s, Math.floor(sReq / 25) * 25));
  }

  const Vs = (Av * fyt * d) / (s * 1000);
  const phiVn = phi * (Vc + Vs);
  return { Vu_kN, Vc, phiVc, VsReq, sMax: sMaxCode, s, phiVn, ok: phiVn >= Vu_kN };
}

function drawSection(svgId, b, h, cc, ds, db, nTop, nBot, stirrupText, zoneLabel) {
  const svg = document.getElementById(svgId);
  const scale = Math.min(145 / b, 175 / h);
  const wPx = b * scale;
  const hPx = h * scale;
  const x0 = (290 - wPx) / 2 - 28;
  const y0 = 28;
  const ccPx = cc * scale;
  const r = Math.max(3, (db / 2) * scale);

  const makeBars = (cnt, top) => {
    let out = "";
    const usableW = Math.max(10, wPx - 2 * ccPx - 2 * r);
    const y = top ? y0 + ccPx + r + 2 : y0 + hPx - ccPx - r - 2;
    for (let i = 0; i < cnt; i++) {
      const x = cnt === 1 ? x0 + wPx / 2 : x0 + ccPx + r + (i * usableW) / (cnt - 1);
      out += `<circle cx="${x}" cy="${y}" r="${r}" fill="#dc2626" stroke="#7f1d1d" stroke-width="0.8"/>`;
    }
    return out;
  };

  svg.innerHTML = `
    <rect x="${x0}" y="${y0}" width="${wPx}" height="${hPx}" fill="#e2e8f0" stroke="#0f172a" stroke-width="2"/>
    <rect x="${x0 + ccPx}" y="${y0 + ccPx}" width="${wPx - 2 * ccPx}" height="${hPx - 2 * ccPx}" rx="4" fill="none" stroke="#2563eb" stroke-width="2"/>
    ${makeBars(nTop, true)}
    ${makeBars(nBot, false)}
    <text x="${x0 + wPx / 2}" y="${y0 + hPx + 18}" text-anchor="middle" font-size="11" font-weight="bold" fill="#0f172a">${b} × ${h} mm</text>
    <text x="${x0 + wPx / 2}" y="${y0 - 10}" text-anchor="middle" font-size="10" fill="#475569">Selimut cc = ${cc} mm</text>
    <text x="${x0 + wPx + 10}" y="${y0 + 18}" font-size="10" font-weight="bold" fill="#dc2626">Atas: ${nTop}D${db}</text>
    <text x="${x0 + wPx + 10}" y="${y0 + hPx - 8}" font-size="10" font-weight="bold" fill="#dc2626">Bawah: ${nBot}D${db}</text>
    <text x="${x0 + wPx + 10}" y="${y0 + hPx / 2 - 4}" font-size="10" font-weight="bold" fill="#2563eb">Sengkang:</text>
    <text x="${x0 + wPx + 10}" y="${y0 + hPx / 2 + 10}" font-size="10" font-weight="bold" fill="#2563eb">${stirrupText}</text>
  `;
}

// Sketsa potongan memanjang dengan kerapatan sengkang proporsional terhadap sTump & sLap
function drawLongitudinal(Ln, b, h, ds, db, nTump, nLap, sTump, sLap) {
  const svg = document.getElementById("svgLongitudinal");
  const x0 = 65, w = 650, y0 = 46, bh = 95;
  const Ln_mm = Ln * 1000;
  const zoneLen_mm = Ln_mm * 0.25; // Zona tumpuan 1/4 Ln kiri & kanan
  const zW = w * 0.25;
  const pxPerMm = w / Ln_mm;

  let stirrups = "";

  // Gambar sengkang zona tumpuan kiri (dari x = 50 mm sampai zoneLen_mm)
  for (let pos = 50; pos <= zoneLen_mm; pos += sTump) {
    const x = x0 + pos * pxPerMm;
    stirrups += `<line x1="${x}" y1="${y0 + 6}" x2="${x}" y2="${y0 + bh - 6}" stroke="#1d4ed8" stroke-width="1.5"/>`;
  }

  // Gambar sengkang zona lapangan (dari zoneLen_mm + sLap sampai Ln_mm - zoneLen_mm - sLap)
  for (let pos = zoneLen_mm + sLap; pos < Ln_mm - zoneLen_mm; pos += sLap) {
    const x = x0 + pos * pxPerMm;
    stirrups += `<line x1="${x}" y1="${y0 + 6}" x2="${x}" y2="${y0 + bh - 6}" stroke="#0ea5e9" stroke-width="1.3"/>`;
  }

  // Gambar sengkang zona tumpuan kanan (dari Ln_mm - zoneLen_mm sampai Ln_mm - 50 mm)
  for (let pos = Ln_mm - zoneLen_mm; pos <= Ln_mm - 50; pos += sTump) {
    const x = x0 + pos * pxPerMm;
    stirrups += `<line x1="${x}" y1="${y0 + 6}" x2="${x}" y2="${y0 + bh - 6}" stroke="#1d4ed8" stroke-width="1.5"/>`;
  }

  svg.innerHTML = `
    <!-- Badan Balok -->
    <rect x="${x0}" y="${y0}" width="${w}" height="${bh}" fill="#f1f5f9" stroke="#0f172a" stroke-width="2"/>
    <!-- Tumpuan Kiri & Kanan -->
    <rect x="${x0 - 18}" y="${y0 - 14}" width="18" height="${bh + 28}" fill="#94a3b8" stroke="#334155"/>
    <rect x="${x0 + w}" y="${y0 - 14}" width="18" height="${bh + 28}" fill="#94a3b8" stroke="#334155"/>
    <!-- Garis Batas Zona Tumpuan & Lapangan -->
    <line x1="${x0 + zW}" y1="${y0 - 8}" x2="${x0 + zW}" y2="${y0 + bh + 8}" stroke="#64748b" stroke-dasharray="4"/>
    <line x1="${x0 + w - zW}" y1="${y0 - 8}" x2="${x0 + w - zW}" y2="${y0 + bh + 8}" stroke="#64748b" stroke-dasharray="4"/>
    <!-- Sengkang Proporsional -->
    ${stirrups}
    <!-- Tulangan Utama Atas & Bawah -->
    <line x1="${x0 - 10}" y1="${y0 + 10}" x2="${x0 + w + 10}" y2="${y0 + 10}" stroke="#dc2626" stroke-width="2.6"/>
    <line x1="${x0 - 10}" y1="${y0 + bh - 10}" x2="${x0 + w + 10}" y2="${y0 + bh - 10}" stroke="#dc2626" stroke-width="2.6"/>
    <!-- Anotasi Tulangan Lentur -->
    <text x="${x0 + zW / 2}" y="${y0 - 14}" text-anchor="middle" font-size="11" font-weight="bold" fill="#dc2626">Atas: ${nTump}D${db}</text>
    <text x="${x0 + w / 2}" y="${y0 - 14}" text-anchor="middle" font-size="11" font-weight="bold" fill="#dc2626">Atas: ${nLap}D${db}</text>
    <text x="${x0 + w - zW / 2}" y="${y0 - 14}" text-anchor="middle" font-size="11" font-weight="bold" fill="#dc2626">Atas: ${nTump}D${db}</text>
    <text x="${x0 + w / 2}" y="${y0 + bh + 16}" text-anchor="middle" font-size="11" font-weight="bold" fill="#dc2626">Tulangan Bawah: ${nLap}D${db} (Menerus)</text>
    <!-- Anotasi Sengkang -->
    <text x="${x0 + zW / 2}" y="${y0 + bh + 34}" text-anchor="middle" font-size="10" font-weight="bold" fill="#1d4ed8">Tumpuan: Ø${ds}-${sTump}</text>
    <text x="${x0 + w / 2}" y="${y0 + bh + 34}" text-anchor="middle" font-size="10" font-weight="bold" fill="#0284c7">Lapangan: Ø${ds}-${sLap}</text>
    <text x="${x0 + w - zW / 2}" y="${y0 + bh + 34}" text-anchor="middle" font-size="10" font-weight="bold" fill="#1d4ed8">Tumpuan: Ø${ds}-${sTump}</text>
    <!-- Keterangan Bentang -->
    <text x="${x0 + w / 2}" y="${y0 + bh + 54}" text-anchor="middle" font-size="11" font-weight="bold" fill="#0f172a">Bentang Ln = ${fmt(Ln, 2)} m (Zona Tumpuan = ${fmt(Ln * 0.25, 2)} m | Zona Lapangan = ${fmt(Ln * 0.5, 2)} m)</text>
  `;
}

function drawCharts(forces) {
  if (sfdChart) sfdChart.destroy();
  if (bmdChart) bmdChart.destroy();

  sfdChart = new Chart(document.getElementById("chartSFD"), {
    type: "line",
    data: {
      labels: forces.x,
      datasets: [{ label: "Geser V (kN)", data: forces.V, borderColor: "#2563eb", fill: true, backgroundColor: "rgba(37,99,235,0.12)", pointRadius: 0 }]
    },
    options: { responsive: true, maintainAspectRatio: false, plugins: { title: { display: true, text: "Diagram Gaya Geser (SFD)" } } }
  });

  bmdChart = new Chart(document.getElementById("chartBMD"), {
    type: "line",
    data: {
      labels: forces.x,
      datasets: [{ label: "Momen M (kN.m)", data: forces.M, borderColor: "#dc2626", fill: true, backgroundColor: "rgba(220,38,38,0.12)", pointRadius: 0 }]
    },
    options: { responsive: true, maintainAspectRatio: false, plugins: { title: { display: true, text: "Diagram Momen Lentur (BMD)" } }, scales: { y: { reverse: true } } }
  });
}

function runCalc() {
  const Ln = parseFloat(document.getElementById("ln").value);
  const support = document.getElementById("support").value;
  const qDL = parseFloat(document.getElementById("qDL").value);
  const qLL = parseFloat(document.getElementById("qLL").value);
  const cc = parseFloat(document.getElementById("cc").value);
  const fc = parseFloat(document.getElementById("fc").value);
  const gammaC = parseFloat(document.getElementById("gammaC").value);
  const fy = parseFloat(document.getElementById("fy").value);
  const db = parseFloat(document.getElementById("db").value);
  const fyt = parseFloat(document.getElementById("fyt").value);
  const ds = parseFloat(document.getElementById("ds").value);

  if (!Ln || !support || Number.isNaN(qDL) || Number.isNaN(qLL) || !cc || !fc || !gammaC || !fy || !db || !fyt || !ds) {
    alert("Lengkapi seluruh data yang diketahui terlebih dahulu sebelum menekan Hitung.");
    return;
  }

  // 1. Estimasi Dimensi Penampang Otomatis
  const dim = estimateDimensions(Ln, support, qDL, qLL, fc, gammaC, fy, cc, ds, db);
  const { b, h, d, hMin, formula } = dim;

  // 2. Analisis Pembebanan & Gaya Dalam
  const qSw = (b / 1000) * (h / 1000) * gammaC;
  const D = qDL + qSw;
  const qu = 1.2 * D + 1.6 * qLL;
  const forces = calcForces(Ln, qu, support);
  const Vu_d = Math.max(0, forces.VuMax - qu * (d / 1000));
  // Gaya geser pada awal zona lapangan (x = 0,25 Ln dari tumpuan)
  const Vu_lap = Math.max(0, forces.VuMax - qu * (0.25 * Ln));

  // 3. Desain Tulangan Lentur
  const flexTump = designFlexure(forces.MuNeg, b, d, fc, fy, db);
  const flexLap = designFlexure(forces.MuPos, b, d, fc, fy, db);

  // 4. Desain Tulangan Geser (Tumpuan vs Lapangan)
  const shearTump = designShear(Vu_d, b, d, fc, fyt, ds, "tumpuan");
  const shearLap = designShear(Vu_lap, b, d, fc, fyt, ds, "lapangan");

  const allOk = h >= hMin && flexTump.ok && flexLap.ok && shearTump.ok && shearLap.ok;
  const badge = (ok) => `<span class="tag ${ok ? "ok" : "bad"}">${ok ? "AMAN" : "TIDAK AMAN"}</span>`;

  document.getElementById("outputSection").classList.remove("hidden");

  const statusBox = document.getElementById("statusBox");
  statusBox.className = `status-bar ${allOk ? "ok" : "bad"}`;
  document.getElementById("statusTitle").textContent = allOk ? "MEMENUHI SNI 2847:2019" : "PERLU PENYESUAIAN";
  document.getElementById("statusText").textContent = allOk
    ? "Estimasi dimensi penampang, tulangan lentur, dan tulangan geser aman."
    : "Terdapat parameter yang belum memenuhi batas SNI.";

  // Output 1: Dimensi
  document.getElementById("outHmin").textContent = `${fmt(hMin, 1)} mm`;
  document.getElementById("outHminFormula").textContent = formula;
  document.getElementById("outDim").textContent = `${b} × ${h} mm`;
  document.getElementById("outEffD").textContent = `${fmt(d, 1)} mm`;

  // Output 2: Pembebanan & Gaya Dalam
  document.getElementById("loadForceBody").innerHTML = `
    <tr><td>Berat Sendiri Balok (q<sub>sw</sub>)</td><td>b &times; h &times; &gamma;<sub>c</sub></td><td><strong>${fmt(qSw)} kN/m</strong></td></tr>
    <tr><td>Total Beban Mati (D)</td><td>q<sub>DL</sub> + q<sub>sw</sub></td><td><strong>${fmt(D)} kN/m</strong></td></tr>
    <tr><td>Beban Terfaktor (q<sub>u</sub>)</td><td>1,2D + 1,6L (SNI 2847:2019 Pasal 5.3.1)</td><td><strong>${fmt(qu)} kN/m</strong></td></tr>
    <tr><td>Momen Tumpuan (M<sub>u</sub><sup>-</sup>)</td><td>Analisis Struktur Tumpuan</td><td><strong>${fmt(forces.MuNeg)} kN.m</strong></td></tr>
    <tr><td>Momen Lapangan (M<sub>u</sub><sup>+</sup>)</td><td>Analisis Struktur Tumpuan</td><td><strong>${fmt(forces.MuPos)} kN.m</strong></td></tr>
    <tr><td>Gaya Geser Maksimum (V<sub>u</sub>) &amp; Kritis (V<sub>u,d</sub>)</td><td>Muka tumpuan &amp; jarak d dari tumpuan</td><td><strong>${fmt(forces.VuMax)} kN</strong> | <strong>${fmt(Vu_d)} kN</strong></td></tr>
  `;

  // Output 3: Tulangan Lentur
  document.getElementById("flexureTableBody").innerHTML = `
    <tr><td>Momen Perlu (M<sub>u</sub>)</td><td>${fmt(forces.MuNeg)} kN.m</td><td>${fmt(forces.MuPos)} kN.m</td></tr>
    <tr><td>Luas Tulangan Perlu (A<sub>s,perlu</sub>)</td><td>${fmt(flexTump.AsReq, 1)} mm²</td><td>${fmt(flexLap.AsReq, 1)} mm²</td></tr>
    <tr><td>Tulangan Terpasang</td><td><strong>${flexTump.n}D${db}</strong> (${fmt(flexTump.AsProv, 1)} mm²)</td><td><strong>${flexLap.n}D${db}</strong> (${fmt(flexLap.AsProv, 1)} mm²)</td></tr>
    <tr><td>Regangan Tarik (&epsilon;<sub>t</sub>) &amp; Faktor &phi;</td><td>${fmt(flexTump.et, 4)} (&phi; = ${fmt(flexTump.phi, 2)})</td><td>${fmt(flexLap.et, 4)} (&phi; = ${fmt(flexLap.phi, 2)})</td></tr>
    <tr><td>Kapasitas Momen Rencana (&phi;M<sub>n</sub>)</td><td><strong>${fmt(flexTump.phiMn)} kN.m</strong></td><td><strong>${fmt(flexLap.phiMn)} kN.m</strong></td></tr>
    <tr><td>Status Pemeriksaan</td><td>${badge(flexTump.ok)}</td><td>${badge(flexLap.ok)}</td></tr>
  `;

  // Output 4: Tulangan Geser
  document.getElementById("shearTableBody").innerHTML = `
    <tr><td>Gaya Geser Perlu (V<sub>u</sub>)</td><td>${fmt(Vu_d)} kN (jarak d)</td><td>${fmt(Vu_lap)} kN (awal lapangan)</td></tr>
    <tr><td>Kapasitas Geser Beton (V<sub>c</sub> &amp; &phi;V<sub>c</sub>)</td><td>${fmt(shearTump.Vc)} kN (&phi;V<sub>c</sub> = ${fmt(shearTump.phiVc)} kN)</td><td>${fmt(shearLap.Vc)} kN (&phi;V<sub>c</sub> = ${fmt(shearLap.phiVc)} kN)</td></tr>
    <tr><td>Kontribusi Sengkang Perlu (V<sub>s,perlu</sub>)</td><td>${fmt(shearTump.VsReq)} kN</td><td>${fmt(shearLap.VsReq)} kN</td></tr>
    <tr><td>Batas Jarak Maksimum SNI (s<sub>max</sub>)</td><td>${fmt(shearTump.sMax, 1)} mm (zona kritis &le; 150 mm)</td><td>${fmt(shearLap.sMax, 1)} mm (d/2)</td></tr>
    <tr><td>Sengkang Terpasang (2 Kaki)</td><td><strong>Ø${ds}-${shearTump.s}</strong></td><td><strong>Ø${ds}-${shearLap.s}</strong></td></tr>
    <tr><td>Kapasitas Geser Rencana (&phi;V<sub>n</sub>)</td><td><strong>${fmt(shearTump.phiVn)} kN</strong></td><td><strong>${fmt(shearLap.phiVn)} kN</strong></td></tr>
    <tr><td>Status Pemeriksaan</td><td>${badge(shearTump.ok)}</td><td>${badge(shearLap.ok)}</td></tr>
  `;

  // Output 5: Sketsa Penulangan & Grafik
  drawSection("svgTumpuan", b, h, cc, ds, db, flexTump.n, flexLap.n, `Ø${ds}-${shearTump.s}`, "Tumpuan");
  drawSection("svgLapangan", b, h, cc, ds, db, flexLap.n, flexLap.n, `Ø${ds}-${shearLap.s}`, "Lapangan");
  drawLongitudinal(Ln, b, h, ds, db, flexTump.n, flexLap.n, shearTump.s, shearLap.s);
  drawCharts(forces);

  lastData = { Ln, b, h, d, fc, fy, fyt, db, ds, qSw, D, qu, forces, Vu_d, flexTump, flexLap, shearTump, shearLap };
}

function resetAll() {
  document.querySelectorAll(".input-col input").forEach((inp) => {
    inp.value = "";
  });
  document.querySelectorAll(".input-col select").forEach((sel) => {
    sel.selectedIndex = 0;
  });
  document.getElementById("outputSection").classList.add("hidden");
  lastData = null;
}

function downloadPDF() {
  if (!lastData || !window.jspdf) {
    alert("Hitung desain terlebih dahulu sebelum mengunduh PDF.");
    return;
  }
  const { jsPDF } = window.jspdf;
  const doc = new jsPDF();
  const r = lastData;
  doc.setFontSize(13);
  doc.text("Laporan Perencanaan Balok Beton Bertulang (SNI 2847:2019)", 14, 18);
  doc.setFontSize(10);
  doc.text(`1. Estimasi Dimensi: ${r.b} x ${r.h} mm (d = ${fmt(r.d, 1)} mm) | Bentang Ln = ${r.Ln} m`, 14, 28);
  doc.text(`2. Pembebanan: q_sw = ${fmt(r.qSw)} kN/m, D = ${fmt(r.D)} kN/m, qu = ${fmt(r.qu)} kN/m`, 14, 35);
  doc.text(`   Gaya Dalam: Mu(-) = ${fmt(r.forces.MuNeg)} kN.m, Mu(+) = ${fmt(r.forces.MuPos)} kN.m, Vu,d = ${fmt(r.Vu_d)} kN`, 14, 42);
  doc.text(`3. Tulangan Lentur: Tumpuan ${r.flexTump.n}D${r.db} (phiMn = ${fmt(r.flexTump.phiMn)} kN.m) | Lapangan ${r.flexLap.n}D${r.db} (phiMn = ${fmt(r.flexLap.phiMn)} kN.m)`, 14, 49);
  doc.text(`4. Tulangan Geser: Tumpuan Ø${r.ds}-${r.shearTump.s} (phiVn = ${fmt(r.shearTump.phiVn)} kN) | Lapangan Ø${r.ds}-${r.shearLap.s} (phiVn = ${fmt(r.shearLap.phiVn)} kN)`, 14, 56);
  doc.save("Laporan_Desain_Balok_SNI2847.pdf");
}

document.addEventListener("DOMContentLoaded", () => {
  document.getElementById("btnCalc").addEventListener("click", runCalc);
  document.getElementById("btnReset").addEventListener("click", resetAll);
  document.getElementById("btnPDF").addEventListener("click", downloadPDF);
  document.getElementById("btnTheme").addEventListener("click", () => {
    const html = document.documentElement;
    html.setAttribute("data-theme", html.getAttribute("data-theme") === "dark" ? "light" : "dark");
  });

  resetAll();
});