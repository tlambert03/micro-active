import { LitElement, html, css } from 'https://cdn.jsdelivr.net/gh/lit/dist@3/core/lit-core.min.js';

// ── Physics utilities ──────────────────────────────────────────────

const PI = Math.PI;

function besselJ1(x) {
  if (x === 0) return 0;
  const ax = Math.abs(x);
  let result;
  if (ax < 8.0) {
    const y = x * x;
    result = x * (72362614232.0 + y * (-7895059235.0 + y *
      (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y *
      (-30.16036606))))));
    result /= (144725228442.0 + y * (2300535178.0 + y *
      (18583304.74 + y * (99447.43394 + y * (376.9991397 + y)))));
  } else {
    const z = 8.0 / ax;
    const y = z * z;
    const xx = ax - 2.356194491;
    const p1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 +
      y * (0.2457520174e-5 + y * (-0.240337019e-6))));
    const q1 = 0.04687499995 + y * (-0.2002690873e-3 +
      y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
    result = Math.sqrt(0.636619772 / ax) *
      (Math.cos(xx) * p1 - z * Math.sin(xx) * q1);
    if (x < 0) result = -result;
  }
  return result;
}

/** Airy disk intensity, normalized so first zero is at r = 1 AU. */
function airyDisk(r) {
  const rho = PI * Math.abs(r) * 1.22;
  if (rho < 1e-12) return 1.0;
  const v = 2 * besselJ1(rho) / rho;
  return v * v;
}

/** Evaluate airyDisk on an array, return new Float64Array. */
function airyDiskArray(grid) {
  const out = new Float64Array(grid.length);
  for (let i = 0; i < grid.length; i++) out[i] = airyDisk(grid[i]);
  return out;
}

/** Top-hat (rect) function: 1 inside radius, 0 outside. */
function tophat(r, radius) {
  return Math.abs(r) < radius ? 1 : 0;
}

/**
 * 1-D convolution of two sampled arrays on a uniform grid.
 * Returns new Float64Array of same length.
 */
function convolve1D(a, b, dx) {
  const n = a.length;
  const out = new Float64Array(n);
  const half = n >> 1;
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let j = 0; j < n; j++) {
      const k = i - j + half;
      if (k >= 0 && k < n) sum += a[j] * b[k];
    }
    out[i] = sum * dx;
  }
  return out;
}

/**
 * Confocal PSF = h_exc(r) × [h_det(r) ⊛ pinhole(r)].
 * Returns { profile, convDet } on the given grid.
 */
function confocalPSF(grid, pinholeRadius) {
  const n = grid.length;
  const dx = grid[1] - grid[0];
  const exc = airyDiskArray(grid);

  if (pinholeRadius > 10) {
    // Wide-open pinhole → widefield: just excitation PSF
    return { profile: exc, convDet: null };
  }

  const det = airyDiskArray(grid);
  const pin = new Float64Array(n);
  for (let i = 0; i < n; i++) pin[i] = tophat(grid[i], pinholeRadius);

  const convDet = convolve1D(det, pin, dx);
  const profile = new Float64Array(n);
  for (let i = 0; i < n; i++) profile[i] = exc[i] * convDet[i];
  return { profile, convDet };
}

/**
 * Off-axis effective PSF: h_exc(r) × h_det(r - offset).
 * 1D version for step 2 (pinhole displaced along the axis).
 */
function offAxisPSF(grid, offset) {
  const n = grid.length;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    out[i] = airyDisk(grid[i]) * airyDisk(grid[i] - offset);
  }
  return out;
}

/**
 * Off-axis effective PSF for a detector element at (dx, dy).
 * Excitation PSF is centered at origin; detection PSF at (dx, dy).
 * We evaluate the product along the 1D x-axis (y=0 slice).
 */
function offAxisPSF2D(grid, dx, dy) {
  const n = grid.length;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const excR = Math.abs(grid[i]); // distance from origin along x-axis
    const detR = Math.sqrt((grid[i] - dx) * (grid[i] - dx) + dy * dy);
    out[i] = airyDisk(excR) * airyDisk(detR);
  }
  return out;
}

/** FWHM of a sampled profile. Returns width in grid units. */
function measureFWHM(profile, grid) {
  let max = 0;
  for (let i = 0; i < profile.length; i++) {
    if (profile[i] > max) max = profile[i];
  }
  const half = max / 2;
  let left = grid[0], right = grid[grid.length - 1];
  for (let i = 1; i < profile.length; i++) {
    if (profile[i - 1] < half && profile[i] >= half) {
      const t = (half - profile[i - 1]) / (profile[i] - profile[i - 1]);
      left = grid[i - 1] + t * (grid[i] - grid[i - 1]);
      break;
    }
  }
  for (let i = profile.length - 2; i >= 0; i--) {
    if (profile[i + 1] < half && profile[i] >= half) {
      const t = (half - profile[i + 1]) / (profile[i] - profile[i + 1]);
      right = grid[i + 1] + t * (grid[i] - grid[i + 1]);
      break;
    }
  }
  return right - left;
}

/** Trapezoidal integration. */
function integrate(profile, grid) {
  let sum = 0;
  for (let i = 1; i < profile.length; i++) {
    sum += (profile[i] + profile[i - 1]) * (grid[i] - grid[i - 1]) / 2;
  }
  return sum;
}

/** Find the position of the peak value. */
function findPeak(profile, grid) {
  let maxVal = -Infinity, maxIdx = 0;
  for (let i = 0; i < profile.length; i++) {
    if (profile[i] > maxVal) { maxVal = profile[i]; maxIdx = i; }
  }
  // Parabolic interpolation around the peak
  if (maxIdx > 0 && maxIdx < profile.length - 1) {
    const a = profile[maxIdx - 1], b = profile[maxIdx], c = profile[maxIdx + 1];
    const denom = 2 * (2 * b - a - c);
    if (Math.abs(denom) > 1e-15) {
      const shift = (a - c) / denom;
      return grid[maxIdx] + shift * (grid[1] - grid[0]);
    }
  }
  return grid[maxIdx];
}

/**
 * Hexagonal detector layout (Airyscan-like, 32 elements).
 * Close-packed hex grid: generate all hex positions within a radius,
 * then keep the 32 closest to center.
 * Returns array of {x, y, ring} in AU.
 */
function hexDetectorLayout() {
  const spacing = 0.35; // AU between adjacent element centers
  const rowH = spacing * Math.sqrt(3) / 2;
  const candidates = [];

  // Generate hex grid candidates over a generous range
  for (let row = -4; row <= 4; row++) {
    const y = row * rowH;
    const xShift = (row % 2 !== 0) ? spacing / 2 : 0;
    for (let col = -4; col <= 4; col++) {
      const x = col * spacing + xShift;
      const dist = Math.sqrt(x * x + y * y);
      candidates.push({ x, y, dist });
    }
  }

  // Sort by distance from center, take closest 32
  candidates.sort((a, b) => a.dist - b.dist);
  return candidates.slice(0, 32).map(({ x, y, dist }) => {
    const ring = dist < spacing * 0.5 ? 0 : dist < spacing * 1.2 ? 1
      : dist < spacing * 2.2 ? 2 : 3;
    return { x, y, ring };
  });
}

/**
 * ISM reconstruction: sum shifted off-axis PSFs for all detector elements.
 * Each element at (dx, dy) contributes a 1D profile shifted back by
 * -|d|/2 (projected onto the 1D axis, using dx component).
 */
/** Signed radial distance: magnitude from (x,y), sign from x. */
function signedRadius(el) {
  const r = Math.sqrt(el.x * el.x + el.y * el.y);
  return el.x < 0 ? -r : r;
}

function ismReconstruction(grid, detectorElements) {
  const n = grid.length;
  const dx = grid[1] - grid[0];
  const result = new Float64Array(n);

  for (const el of detectorElements) {
    const psf = offAxisPSF2D(grid, el.x, el.y);
    const shiftBins = (-el.x / 2) / dx; // reassign by x-projection

    for (let i = 0; i < n; i++) {
      const srcIdx = i - shiftBins;
      const lo = Math.floor(srcIdx);
      const frac = srcIdx - lo;
      if (lo >= 0 && lo + 1 < n) {
        result[i] += psf[lo] * (1 - frac) + psf[lo + 1] * frac;
      }
    }
  }
  return result;
}

// ── Drawing utilities ──────────────────────────────────────────────

/**
 * Draw a labeled 1-D plot on canvas.
 * rect: {x, y, w, h} — pixel region for the plot area (inside axes).
 * curves: [{data: Float64Array, color, lineWidth, dash?, label?}, ...]
 * options: {grid: positions, xRange:[min,max], yRange:[min,max],
 *           xLabel, yLabel, title, annotations?:[{type,…}]}
 */
function drawPlot(ctx, rect, curves, options) {
  const { x, y, w, h } = rect;
  const { grid, xRange, yRange, xLabel, yLabel, title } = options;
  const [xMin, xMax] = xRange;
  const [yMin, yMax] = yRange;

  const toX = (v) => x + (v - xMin) / (xMax - xMin) * w;
  const toY = (v) => y + h - (v - yMin) / (yMax - yMin) * h;

  // Compute nice y-tick interval
  const ySpan = yMax - yMin;
  const rawStep = ySpan / 5;
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const residual = rawStep / mag;
  const yStep = mag * (residual > 5 ? 10 : residual > 2 ? 5 : residual > 1 ? 2 : 1);
  const yDecimals = yStep < 1 ? Math.max(0, -Math.floor(Math.log10(yStep))) + 1 : 0;

  // Grid lines
  ctx.strokeStyle = 'rgba(255,255,255,0.08)';
  ctx.lineWidth = 0.5;
  for (let gx = Math.ceil(xMin); gx <= xMax; gx++) {
    ctx.beginPath();
    ctx.moveTo(toX(gx), y);
    ctx.lineTo(toX(gx), y + h);
    ctx.stroke();
  }
  for (let gy = yStep; gy <= yMax; gy += yStep) {
    ctx.beginPath();
    ctx.moveTo(x, toY(gy));
    ctx.lineTo(x + w, toY(gy));
    ctx.stroke();
  }

  // Axes
  ctx.strokeStyle = 'rgba(255,255,255,0.4)';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(x, y + h);
  ctx.lineTo(x + w, y + h);
  ctx.moveTo(x, y);
  ctx.lineTo(x, y + h);
  ctx.stroke();

  // Tick labels
  ctx.fillStyle = 'rgba(255,255,255,0.5)';
  ctx.font = '11px system-ui';
  ctx.textAlign = 'center';
  for (let gx = Math.ceil(xMin); gx <= xMax; gx++) {
    ctx.fillText(gx.toString(), toX(gx), y + h + 14);
  }
  ctx.textAlign = 'right';
  for (let gy = 0; gy <= yMax + yStep * 0.01; gy += yStep) {
    ctx.fillText(gy.toFixed(yDecimals), x - 4, toY(gy) + 4);
  }

  // Axis labels
  ctx.fillStyle = 'rgba(255,255,255,0.6)';
  ctx.font = '12px system-ui';
  ctx.textAlign = 'center';
  if (xLabel) ctx.fillText(xLabel, x + w / 2, y + h + 30);
  if (yLabel) {
    ctx.save();
    ctx.translate(x - 30, y + h / 2);
    ctx.rotate(-PI / 2);
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();
  }

  // Title
  if (title) {
    ctx.fillStyle = 'rgba(255,255,255,0.8)';
    ctx.font = 'bold 13px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText(title, x + w / 2, y - 8);
  }

  // Curves
  for (const curve of curves) {
    if (curve.legendOnly) continue;
    ctx.strokeStyle = curve.color;
    ctx.lineWidth = curve.lineWidth || 2;
    ctx.setLineDash(curve.dash || []);
    ctx.beginPath();
    let started = false;
    for (let i = 0; i < grid.length; i++) {
      if (curve.xEnd !== undefined && grid[i] > curve.xEnd) break;
      const px = toX(grid[i]);
      const py = toY(curve.data[i]);
      if (px >= x && px <= x + w) {
        if (!started) { ctx.moveTo(px, py); started = true; }
        else ctx.lineTo(px, py);
      }
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Legend
  const legendCurves = curves.filter(c => c.label);
  if (legendCurves.length) {
    const lx = x + w - 10;
    let ly = y + 16;
    ctx.textAlign = 'right';
    ctx.font = '11px system-ui';
    for (const c of legendCurves) {
      ctx.strokeStyle = c.color;
      ctx.lineWidth = c.lineWidth || 2;
      ctx.setLineDash(c.dash || []);
      ctx.beginPath();
      ctx.moveTo(lx - 40, ly);
      ctx.lineTo(lx - 22, ly);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = c.color;
      ctx.fillText(c.label, lx - 44, ly + 4);
      ly += 16;
    }
  }

  // Annotations
  if (options.annotations) {
    for (const ann of options.annotations) {
      if (ann.type === 'vline') {
        ctx.strokeStyle = ann.color || 'rgba(255,255,255,0.4)';
        ctx.lineWidth = 1;
        ctx.setLineDash(ann.dash || [4, 4]);
        const px = toX(ann.x);
        ctx.beginPath();
        ctx.moveTo(px, y);
        ctx.lineTo(px, y + h);
        ctx.stroke();
        ctx.setLineDash([]);
        if (ann.label) {
          ctx.fillStyle = ann.color || 'rgba(255,255,255,0.6)';
          ctx.font = '11px system-ui';
          ctx.textAlign = 'center';
          ctx.fillText(ann.label, px, y + h - 6);
        }
      }
      if (ann.type === 'fwhm') {
        const halfMax = ann.halfMax;
        const py = toY(halfMax);
        ctx.strokeStyle = ann.color || '#fff';
        ctx.lineWidth = 1;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(toX(ann.left), py);
        ctx.lineTo(toX(ann.right), py);
        ctx.stroke();
        ctx.setLineDash([]);
        // Arrow heads
        for (const xp of [ann.left, ann.right]) {
          ctx.beginPath();
          ctx.arc(toX(xp), py, 3, 0, 2 * PI);
          ctx.fillStyle = ann.color || '#fff';
          ctx.fill();
        }
        if (ann.label) {
          ctx.fillStyle = ann.color || '#fff';
          ctx.font = '11px system-ui';
          ctx.textAlign = 'center';
          ctx.fillText(ann.label, toX((ann.left + ann.right) / 2), py - 6);
        }
      }
    }
  }

  return { toX, toY };
}

// ── Component ──────────────────────────────────────────────────────

class IsmTutorial extends LitElement {
  static properties = {
    step: { type: Number },
    pinholeSize: { type: Number, attribute: 'pinhole-size' },
    normalized: { type: Boolean },
    pinholeOffset: { type: Number, attribute: 'pinhole-offset' },
    pinholeSize2: { type: Number, attribute: 'pinhole-size-2' },
    beamPos: { type: Number, attribute: 'beam-pos' },
    pinholeSize4: { type: Number, attribute: 'pinhole-size-4' },
    pinholeOffset4: { type: Number, attribute: 'pinhole-offset-4' },
    stokesShift4: { type: Number, attribute: 'stokes-shift-4' },
    selectedElement: { type: Number, attribute: 'selected-element' },
    showReassignment: { type: Boolean, attribute: 'show-reassignment' },
    _reassignProgress: { state: true },
    _showConvDet: { state: true },
    _shifted: { state: true },
  };

  static styles = css`
    :host {
      display: block;
      font-family: system-ui, -apple-system, sans-serif;
      max-width: 900px;
    }

    .container {
      background: #111122;
      border-radius: 10px;
      overflow: hidden;
    }

    .step-nav {
      display: flex;
      gap: 2px;
      background: #0a0a18;
      padding: 2px;
    }

    .step-btn {
      flex: 1;
      padding: 10px 8px;
      background: transparent;
      border: none;
      color: rgba(255,255,255,0.4);
      font: 600 13px system-ui;
      cursor: pointer;
      transition: all 0.2s;
      text-align: center;
    }

    .step-btn:hover {
      color: rgba(255,255,255,0.7);
      background: rgba(255,255,255,0.04);
    }

    .step-btn.active {
      color: #7eb8ff;
      background: rgba(126,184,255,0.1);
      border-bottom: 2px solid #7eb8ff;
    }

    .step-description {
      padding: 12px 16px;
      color: rgba(255,255,255,0.75);
      font-size: 13.5px;
      line-height: 1.6;
      background: rgba(255,255,255,0.03);
      border-bottom: 1px solid rgba(255,255,255,0.06);
    }

    .step-description strong {
      color: #7eb8ff;
    }

    canvas {
      display: block;
      width: 100%;
      background: #0d0d20;
    }

    .controls {
      padding: 14px 16px;
      display: flex;
      flex-wrap: wrap;
      gap: 12px 20px;
      align-items: center;
      background: rgba(255,255,255,0.03);
      border-top: 1px solid rgba(255,255,255,0.06);
    }

    .control-group {
      display: flex;
      align-items: center;
      gap: 8px;
      min-width: 200px;
      flex: 1;
    }

    .control-group label {
      color: rgba(255,255,255,0.6);
      font-size: 12px;
      font-weight: 500;
      white-space: nowrap;
    }

    sl-range {
      flex: 1;
      --track-color-active: #4a7ab5;
      --track-color-inactive: #222238;
      --track-height: 4px;
      --thumb-size: 18px;
      --sl-color-primary-600: #7eb8ff;
      --sl-color-primary-700: #6aa8ee;
      --sl-input-background-color-hover: #7eb8ff;
      --tooltip-offset: -9999px;
      min-width: 120px;
      cursor: pointer;
    }

    sl-range::part(tooltip) {
      display: none;
    }

    sl-range::part(input)::-webkit-slider-thumb {
      background: #7eb8ff;
    }

    sl-range::part(input)::-webkit-slider-thumb:hover {
      background: #9ecbff;
    }

    .value {
      color: #7eb8ff;
      font: 600 12px monospace;
      min-width: 52px;
      text-align: right;
    }

    .btn-group {
      display: flex;
      gap: 2px;
    }

    .btn-group button {
      padding: 6px 12px;
      background: rgba(255,255,255,0.06);
      border: 1px solid rgba(255,255,255,0.1);
      color: rgba(255,255,255,0.6);
      font: 500 12px system-ui;
      cursor: pointer;
      transition: all 0.15s;
    }

    .btn-group button:first-child { border-radius: 4px 0 0 4px; }
    .btn-group button:last-child { border-radius: 0 4px 4px 0; }
    .btn-group button:only-child { border-radius: 4px; }

    .btn-group button.active {
      background: rgba(126,184,255,0.2);
      color: #7eb8ff;
      border-color: rgba(126,184,255,0.3);
    }

    .btn-group button:hover:not(.active) {
      background: rgba(255,255,255,0.1);
      color: rgba(255,255,255,0.8);
    }

    .info-bar {
      padding: 8px 16px;
      color: rgba(255,255,255,0.5);
      font: 12px monospace;
      background: rgba(0,0,0,0.2);
      display: flex;
      gap: 20px;
      flex-wrap: wrap;
    }

    .info-bar span {
      white-space: nowrap;
    }

    .info-bar .highlight {
      color: #7eb8ff;
    }
  `;

  constructor() {
    super();
    this.step = 1;
    this.pinholeSize = 1.0;
    this.normalized = true;
    this.pinholeOffset = 0.8;
    this.pinholeSize2 = 0;
    this.beamPos = -0.4;
    this.pinholeSize4 = 0.5;
    this.pinholeOffset4 = 0;
    this.stokesShift4 = 60;
    this.selectedElement = -1; // -1 = all
    this._lastSingleElement = 0;
    this.showReassignment = false;
    this._reassignProgress = 1.0;
    this._showConvDet = true;
    this._shifted = false;

    // Pre-computed data
    this._grid = null;
    this._fwhmCurve = null;
    this._signalCurve = null;
    this._pinholeSweep = null;
    this._detectorElements = hexDetectorLayout();
    this._animFrameId = null;
  }

  connectedCallback() {
    super.connectedCallback();
    this._initGrid();
    this._precompute();
  }

  _initGrid() {
    const n = 601;
    const range = 3;
    this._grid = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      this._grid[i] = -range + (2 * range * i) / (n - 1);
    }
  }

  _precompute() {
    const grid = this._grid;
    const nSweep = 80;
    this._pinholeSweep = new Float64Array(nSweep);
    this._fwhmCurve = new Float64Array(nSweep);
    this._signalCurve = new Float64Array(nSweep);

    // Widefield reference
    const wfPSF = airyDiskArray(grid);
    const wfFWHM = measureFWHM(wfPSF, grid);
    const wfSignal = integrate(wfPSF, grid);

    for (let i = 0; i < nSweep; i++) {
      const ph = 0.02 + (3.0 - 0.02) * i / (nSweep - 1);
      this._pinholeSweep[i] = ph;
      const { profile } = confocalPSF(grid, ph);
      this._fwhmCurve[i] = measureFWHM(profile, grid) / wfFWHM;
      this._signalCurve[i] = integrate(profile, grid) / wfSignal;
    }
  }

  async firstUpdated() {
    // Load Shoelace
    if (!document.querySelector('link[href*="shoelace"]')) {
      const link = document.createElement('link');
      link.rel = 'stylesheet';
      link.href = 'https://cdn.jsdelivr.net/npm/@shoelace-style/shoelace@2.12.0/cdn/themes/dark.css';
      document.head.appendChild(link);
    }
    if (!customElements.get('sl-range')) {
      await import('https://cdn.jsdelivr.net/npm/@shoelace-style/shoelace@2.12.0/cdn/components/range/range.js');
    }
    this._canvas = this.shadowRoot.querySelector('canvas');
    this._ctx = this._canvas.getContext('2d');
    this._renderFrame();
  }

  updated(changed) {
    if (this._ctx) {
      if (changed.has('showReassignment') && this.showReassignment) {
        this._startReassignmentAnimation();
      } else {
        this._renderFrame();
      }
    }
  }

  // ── Rendering dispatch ─────────────────────────────────────────

  _renderFrame() {
    const ctx = this._ctx;
    if (!ctx) return;
    const canvas = this._canvas;
    // Handle DPR for crisp rendering
    const dpr = window.devicePixelRatio || 1;
    const cssW = canvas.clientWidth || 860;
    const cssH = 460;
    canvas.width = cssW * dpr;
    canvas.height = cssH * dpr;
    canvas.style.height = cssH + 'px';
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    ctx.clearRect(0, 0, cssW, cssH);

    if (this.step === 1) this._renderStep1(ctx, cssW, cssH);
    else if (this.step === 2) this._renderStep2(ctx, cssW, cssH);
    else if (this.step === 3) this._renderStep3(ctx, cssW, cssH);
    else if (this.step === 4) this._renderStep4(ctx, cssW, cssH);
  }

  // ── Step 1: Confocal pinhole ───────────────────────────────────

  _renderStep1(ctx, W, H) {
    const grid = this._grid;
    const { profile, convDet } = confocalPSF(grid, this.pinholeSize);
    const exc = airyDiskArray(grid);

    // Normalize or not
    let yMax = 1.1;
    let plotProfile = profile;
    let plotExc = exc;
    let plotDet = convDet;

    if (this.normalized) {
      // Normalize to peak = 1
      let pMax = 0;
      for (let i = 0; i < profile.length; i++) {
        if (profile[i] > pMax) pMax = profile[i];
      }
      if (pMax > 0) {
        plotProfile = new Float64Array(profile.length);
        for (let i = 0; i < profile.length; i++) plotProfile[i] = profile[i] / pMax;
      }
    } else {
      // Unnormalized — scale relative to widefield peak
      let eMax = 0;
      for (let i = 0; i < exc.length; i++) {
        if (exc[i] > eMax) eMax = exc[i];
      }
      // Profile is already in absolute units; find max for y-axis
      let pMax = 0;
      for (let i = 0; i < profile.length; i++) {
        if (profile[i] > pMax) pMax = profile[i];
      }
      yMax = Math.max(1.1, eMax * 1.1);
    }

    // Main plot area
    const mainRect = { x: 50, y: 30, w: W * 0.58, h: H - 70 };

    // Pinhole top-hat
    const pinholeData = new Float64Array(grid.length);
    for (let i = 0; i < grid.length; i++) {
      pinholeData[i] = Math.abs(grid[i]) < this.pinholeSize ? 1 : 0;
    }

    const curves = [
      { data: plotExc, color: '#5599ff', lineWidth: 1.5, dash: [6, 4], label: 'Excitation PSF' },
      { data: pinholeData, color: '#ff776644', lineWidth: 1.5, label: 'Pinhole' },
    ];

    if (this._showConvDet && convDet) {
      let dMax = 0;
      for (let i = 0; i < convDet.length; i++) {
        if (convDet[i] > dMax) dMax = convDet[i];
      }
      const normDet = new Float64Array(convDet.length);
      if (this.normalized && dMax > 0) {
        for (let i = 0; i < convDet.length; i++) normDet[i] = convDet[i] / dMax;
      } else {
        for (let i = 0; i < convDet.length; i++) normDet[i] = convDet[i];
      }
      curves.push({
        data: normDet, color: '#ff7766', lineWidth: 1.5,
        dash: [6, 4], label: 'Detection PSF'
      });
    }

    curves.push(
      { data: plotProfile, color: '#66ff99', lineWidth: 2.5, label: 'Confocal PSF' },
    );

    // FWHM annotation
    const fwhm = measureFWHM(plotProfile, grid);
    let halfMax = 0;
    for (const v of plotProfile) { if (v > halfMax) halfMax = v; }
    halfMax /= 2;

    const annotations = [
      { type: 'fwhm', left: -fwhm / 2, right: fwhm / 2, halfMax,
        color: 'rgba(102,255,153,0.7)',
        label: `FWHM: ${fwhm.toFixed(3)} AU` },
    ];

    drawPlot(ctx, mainRect, curves, {
      grid, xRange: [-3, 3], yRange: [0, yMax],
      xLabel: 'Position (AU)', yLabel: 'Intensity',
      title: 'Confocal PSF vs Pinhole Size',
      annotations,
    });

    // ── Inset: FWHM & Signal vs pinhole ──
    this._drawInsetCurves(ctx, W, H);
  }

  _drawInsetCurves(ctx, W, H) {
    const sweep = this._pinholeSweep;
    const fwhmC = this._fwhmCurve;
    const sigC = this._signalCurve;
    const n = sweep.length;

    const ix = W * 0.65;
    const iy = 40;
    const iw = W * 0.30;
    const ih = (H - 90) / 2 - 10;

    // FWHM inset
    this._drawMiniPlot(ctx, ix, iy, iw, ih,
      sweep, fwhmC, '#66ff99', 'Relative FWHM',
      [0, 3], [0.5, 1.1], this.pinholeSize);

    // Add √2 reference line
    const y707 = iy + ih - (1 / Math.sqrt(2) - 0.5) / (1.1 - 0.5) * ih;
    ctx.strokeStyle = 'rgba(255,255,255,0.2)';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(ix, y707);
    ctx.lineTo(ix + iw, y707);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(255,255,255,0.35)';
    ctx.font = '10px system-ui';
    ctx.textAlign = 'left';
    ctx.fillText('1/√2 ≈ 0.707', ix + 2, y707 - 3);

    // Signal inset
    const sy = iy + ih + 30;
    this._drawMiniPlot(ctx, ix, sy, iw, ih,
      sweep, sigC, '#ffaa55', 'Relative Signal',
      [0, 3], [0, 1.1], this.pinholeSize);
  }

  _drawMiniPlot(ctx, x, y, w, h, xData, yData, color, title, xRange, yRange, markerX) {
    const [xMin, xMax] = xRange;
    const [yMin, yMax] = yRange;
    const toX = (v) => x + (v - xMin) / (xMax - xMin) * w;
    const toY = (v) => y + h - (v - yMin) / (yMax - yMin) * h;

    // Background
    ctx.fillStyle = 'rgba(0,0,0,0.3)';
    ctx.strokeStyle = 'rgba(255,255,255,0.1)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.roundRect(x - 4, y - 18, w + 8, h + 28, 4);
    ctx.fill();
    ctx.stroke();

    // Title
    ctx.fillStyle = 'rgba(255,255,255,0.5)';
    ctx.font = '11px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText(title, x + w / 2, y - 5);

    // Axes
    ctx.strokeStyle = 'rgba(255,255,255,0.2)';
    ctx.beginPath();
    ctx.moveTo(x, y + h);
    ctx.lineTo(x + w, y + h);
    ctx.moveTo(x, y);
    ctx.lineTo(x, y + h);
    ctx.stroke();

    // Tick labels
    ctx.fillStyle = 'rgba(255,255,255,0.35)';
    ctx.font = '9px system-ui';
    ctx.textAlign = 'center';
    for (let t = 0; t <= 3; t++) ctx.fillText(t, toX(t), y + h + 10);
    ctx.textAlign = 'right';
    for (let t = yMin; t <= yMax; t += 0.5) {
      if (t >= yMin) ctx.fillText(t.toFixed(1), x - 2, toY(t) + 3);
    }

    // Curve
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let i = 0; i < xData.length; i++) {
      const px = toX(xData[i]);
      const py = toY(yData[i]);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Marker dot at current pinhole
    if (markerX !== undefined) {
      // Interpolate y value
      let my = 0;
      for (let i = 1; i < xData.length; i++) {
        if (xData[i] >= markerX) {
          const t = (markerX - xData[i - 1]) / (xData[i] - xData[i - 1]);
          my = yData[i - 1] + t * (yData[i] - yData[i - 1]);
          break;
        }
      }
      ctx.beginPath();
      ctx.arc(toX(markerX), toY(my), 5, 0, 2 * PI);
      ctx.fillStyle = color;
      ctx.fill();
      ctx.strokeStyle = '#fff';
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }
  }

  // ── Step 2: Off-axis pinhole ───────────────────────────────────

  _renderStep2(ctx, W, H) {
    const grid = this._grid;
    const offset = this.pinholeOffset;

    const exc = airyDiskArray(grid);
    const n = grid.length;
    const dx = grid[1] - grid[0];
    const phSize = this.pinholeSize2;

    // Detection PSF: Airy disk centered at offset, convolved with pinhole if finite
    const det = new Float64Array(n);
    let eff;
    if (phSize > 0.01) {
      // Finite pinhole: convolve detection Airy with top-hat, then shift
      const detRaw = airyDiskArray(grid);
      const pin = new Float64Array(n);
      for (let i = 0; i < n; i++) pin[i] = tophat(grid[i], phSize);
      const convDet = convolve1D(detRaw, pin, dx);
      // Shift convDet to be centered at offset
      for (let i = 0; i < n; i++) {
        const srcIdx = (grid[i] - offset - grid[0]) / dx;
        const lo = Math.floor(srcIdx);
        const frac = srcIdx - lo;
        det[i] = (lo >= 0 && lo + 1 < n)
          ? convDet[lo] * (1 - frac) + convDet[lo + 1] * frac : 0;
      }
      eff = new Float64Array(n);
      for (let i = 0; i < n; i++) eff[i] = exc[i] * det[i];
    } else {
      // Infinitesimal pinhole: pointwise product
      for (let i = 0; i < n; i++) det[i] = airyDisk(grid[i] - offset);
      eff = offAxisPSF(grid, offset);
    }

    const fwhm = measureFWHM(eff, grid);

    // Reference FWHM (on-axis, same pinhole size)
    let refEff;
    if (phSize > 0.01) {
      const detRaw = airyDiskArray(grid);
      const pin = new Float64Array(n);
      for (let i = 0; i < n; i++) pin[i] = tophat(grid[i], phSize);
      const convDet = convolve1D(detRaw, pin, dx);
      refEff = new Float64Array(n);
      for (let i = 0; i < n; i++) refEff[i] = exc[i] * convDet[i];
    } else {
      refEff = offAxisPSF(grid, 0);
    }
    const refFWHM = measureFWHM(refEff, grid);

    let plotExc = exc;
    let plotDet = det;
    let plotEff = eff;
    let yMax = 1.1;
    let yLabel = 'Intensity';

    if (this.normalized) {
      let eMax = 0;
      for (const v of eff) { if (v > eMax) eMax = v; }
      if (eMax > 0) {
        plotEff = new Float64Array(eff.length);
        for (let i = 0; i < eff.length; i++) plotEff[i] = eff[i] / eMax;
      }
      yLabel = 'Intensity (norm.)';
    } else {
      // Unnormalized: all curves in absolute units (excitation peak = 1)
      yLabel = 'Intensity';
    }

    const mainRect = { x: 50, y: 30, w: W - 100, h: H * 0.62 };

    const curves = [
      { data: plotExc, color: '#5599ff', lineWidth: 1.5, dash: [6, 4], label: 'Excitation PSF' },
      { data: plotDet, color: '#ff7766', lineWidth: 1.5, dash: [6, 4], label: 'Detection PSF' },
      { data: plotEff, color: '#66ff99', lineWidth: 2.5, label: 'Effective PSF' },
    ];

    const reassignPos = offset / 2;

    const annotations = [
      { type: 'vline', x: 0, color: 'rgba(85,153,255,0.4)', dash: [3, 3], label: 'Excitation' },
      { type: 'vline', x: offset, color: 'rgba(255,119,102,0.4)', dash: [3, 3], label: 'Pinhole' },
      { type: 'vline', x: reassignPos, color: 'rgba(102,255,153,0.5)', dash: [3, 3], label: 'Reassign' },
    ];

    drawPlot(ctx, mainRect, curves, {
      grid, xRange: [-3, 3], yRange: [0, yMax],
      xLabel: 'Position (AU)', yLabel,
      title: 'Off-Axis Pinhole Effect',
      annotations,
    });

    // Draw the half-shift annotation below the plot
    const schY = mainRect.y + mainRect.h + 30;
    const schH = H - schY - 10;
    const toX = (v) => mainRect.x + (v - (-3)) / 6 * mainRect.w;

    // Optical axis line
    ctx.strokeStyle = 'rgba(255,255,255,0.2)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(mainRect.x, schY + schH / 2);
    ctx.lineTo(mainRect.x + mainRect.w, schY + schH / 2);
    ctx.stroke();

    // Excitation marker
    ctx.fillStyle = '#5599ff';
    ctx.beginPath();
    ctx.arc(toX(0), schY + schH / 2, 6, 0, 2 * PI);
    ctx.fill();
    ctx.fillStyle = 'rgba(255,255,255,0.6)';
    ctx.font = '11px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText('Excitation', toX(0), schY + schH / 2 - 12);

    // Pinhole marker
    ctx.fillStyle = '#ff7766';
    ctx.beginPath();
    ctx.arc(toX(offset), schY + schH / 2, 6, 0, 2 * PI);
    ctx.fill();
    ctx.fillText('Pinhole', toX(offset), schY + schH / 2 - 12);

    // Reassignment position marker (always exactly d/2)
    ctx.fillStyle = '#66ff99';
    ctx.beginPath();
    ctx.arc(toX(reassignPos), schY + schH / 2, 6, 0, 2 * PI);
    ctx.fill();
    ctx.fillText('Reassign', toX(reassignPos), schY + schH / 2 - 12);

    // Distance annotations
    if (Math.abs(offset) > 0.05) {
      // Arrow from excitation to pinhole
      this._drawDistArrow(ctx, toX(0), toX(offset), schY + schH / 2 + 18,
        `d = ${Math.abs(offset).toFixed(2)} AU`, 'rgba(255,255,255,0.5)');
      // Arrow from excitation to reassignment position
      this._drawDistArrow(ctx, toX(0), toX(reassignPos), schY + schH / 2 + 36,
        `d/2 = ${Math.abs(reassignPos).toFixed(2)} AU`, '#66ff99');
    }

    // Info: FWHM comparison
    ctx.fillStyle = 'rgba(255,255,255,0.4)';
    ctx.font = '11px system-ui';
    ctx.textAlign = 'right';
    ctx.fillText(
      `FWHM: ${fwhm.toFixed(3)} AU  (on-axis ref: ${refFWHM.toFixed(3)} AU)`,
      mainRect.x + mainRect.w, schY + schH - 2
    );
  }

  _drawDistArrow(ctx, x1, x2, y, label, color) {
    ctx.strokeStyle = color;
    ctx.fillStyle = color;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(x1, y);
    ctx.lineTo(x2, y);
    ctx.stroke();
    // Arrow heads
    const dir = x2 > x1 ? 1 : -1;
    ctx.beginPath();
    ctx.moveTo(x2, y);
    ctx.lineTo(x2 - dir * 6, y - 3);
    ctx.lineTo(x2 - dir * 6, y + 3);
    ctx.fill();
    ctx.beginPath();
    ctx.moveTo(x1, y);
    ctx.lineTo(x1 + dir * 6, y - 3);
    ctx.lineTo(x1 + dir * 6, y + 3);
    ctx.fill();
    // Label
    ctx.font = '11px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText(label, (x1 + x2) / 2, y - 5);
  }

  // ── Step 3: ISM reconstruction ─────────────────────────────────

  _renderStep3(ctx, W, H) {
    const grid = this._grid;
    const elements = this._detectorElements;

    // Draw hex detector on the left
    const hexCX = 130;
    const hexCY = H / 2 - 10;
    const hexScale = 95; // pixels per AU
    this._drawHexDetector(ctx, hexCX, hexCY, hexScale);

    // Draw PSF comparison on the right
    const plotRect = { x: 290, y: 30, w: W - 330, h: H - 80 };
    this._drawStep3Plot(ctx, plotRect);
  }

  _drawHexDetector(ctx, cx, cy, scale) {
    const elements = this._detectorElements;
    const radius = scale * 0.18;

    // Title
    ctx.fillStyle = 'rgba(255,255,255,0.6)';
    ctx.font = 'bold 12px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText('Detector Array', cx, 20);
    ctx.font = '11px system-ui';
    ctx.fillStyle = 'rgba(255,255,255,0.35)';
    ctx.fillText(`(${elements.length} elements)`, cx, 34);

    for (let i = 0; i < elements.length; i++) {
      const el = elements[i];
      const ex = cx + el.x * scale;
      const ey = cy + el.y * scale;

      const isSelected = this.selectedElement === i;
      const isAll = this.selectedElement === -1;

      // Hexagon path
      ctx.beginPath();
      for (let v = 0; v < 6; v++) {
        const a = PI / 3 * v - PI / 6;
        const hx = ex + radius * Math.cos(a);
        const hy = ey + radius * Math.sin(a);
        if (v === 0) ctx.moveTo(hx, hy);
        else ctx.lineTo(hx, hy);
      }
      ctx.closePath();

      // Fill color based on selection
      if (isSelected) {
        ctx.fillStyle = '#7eb8ff';
      } else if (isAll) {
        // Color by signal strength (brighter = closer to center)
        const dist = Math.sqrt(el.x * el.x + el.y * el.y);
        const brightness = Math.max(0.15, 1 - dist / 1.5);
        ctx.fillStyle = `rgba(126,184,255,${brightness * 0.5})`;
      } else {
        ctx.fillStyle = 'rgba(255,255,255,0.08)';
      }
      ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,0.25)';
      ctx.lineWidth = 1;
      ctx.stroke();

      // Element index
      if (radius > 12) {
        ctx.fillStyle = isSelected ? '#000' : 'rgba(255,255,255,0.35)';
        ctx.font = '8px system-ui';
        ctx.textAlign = 'center';
        ctx.fillText(i, ex, ey + 3);
      }
    }

    // "All" button below
    const allY = cy + scale * 1.5 + 10;
    ctx.fillStyle = this.selectedElement === -1
      ? 'rgba(126,184,255,0.2)' : 'rgba(255,255,255,0.06)';
    ctx.strokeStyle = this.selectedElement === -1
      ? 'rgba(126,184,255,0.4)' : 'rgba(255,255,255,0.15)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.roundRect(cx - 30, allY, 60, 24, 4);
    ctx.fill();
    ctx.stroke();
    ctx.fillStyle = this.selectedElement === -1 ? '#7eb8ff' : 'rgba(255,255,255,0.5)';
    ctx.font = '11px system-ui';
    ctx.textAlign = 'center';
    ctx.fillText('All (ISM)', cx, allY + 16);

    // Store hit-test data
    this._hexHitData = { cx, cy, scale, radius, elements, allBtnY: allY };
  }

  _drawStep3Plot(ctx, rect) {
    const grid = this._grid;
    const elements = this._detectorElements;
    const norm = this.normalized;

    const wf = airyDiskArray(grid);

    const normArr = (arr) => {
      let m = 0;
      for (const v of arr) { if (v > m) m = v; }
      const out = new Float64Array(arr.length);
      if (m > 0) for (let i = 0; i < arr.length; i++) out[i] = arr[i] / m;
      return out;
    };
    const maybeNorm = (arr) => norm ? normArr(arr) : arr;

    const curves = [];
    let yMax = 1.1;

    const n = grid.length;

    if (this.selectedElement === -1) {
      const dx = grid[1] - grid[0];
      const shifted = this._shifted;

      // Compute sum and individual element curves using radial distance
      const sum = new Float64Array(n);
      const elementCurves = [];
      for (const el of elements) {
        const psf = offAxisPSF2D(grid, el.x, el.y);
        const displayed = new Float64Array(n);

        if (shifted) {
          // Pixel reassignment: shift back by -x/2
          const shiftBins = (-el.x / 2) / dx;
          for (let i = 0; i < n; i++) {
            const srcIdx = i - shiftBins;
            const lo = Math.floor(srcIdx);
            const frac = srcIdx - lo;
            if (lo >= 0 && lo + 1 < n) {
              displayed[i] = psf[lo] * (1 - frac) + psf[lo + 1] * frac;
            }
          }
        } else {
          for (let i = 0; i < n; i++) displayed[i] = psf[i];
        }

        for (let i = 0; i < n; i++) sum[i] += displayed[i];
        elementCurves.push(displayed);
      }

      // Normalize sum to peak = 1 for display scaling
      let sumMax = 0;
      for (const v of sum) { if (v > sumMax) sumMax = v; }
      const sumScale = norm && sumMax > 0 ? 1 / sumMax : 1;

      // Find max of individual element curves to scale them visibly
      let elMax = 0;
      for (const ec of elementCurves) {
        for (const v of ec) { if (v > elMax) elMax = v; }
      }
      // Scale individual curves so the tallest reaches ~0.35 of plot height
      const elScale = elMax > 0 ? (0.35 * (norm ? 1 : sumMax)) / elMax : 1;

      // Individual element curves (opacity scales with peak intensity)
      for (const ec of elementCurves) {
        let peak = 0;
        for (const v of ec) { if (v > peak) peak = v; }
        const alpha = elMax > 0 ? 0.08 + 0.25 * (peak / elMax) : 0.15;
        const scaled = new Float64Array(n);
        for (let i = 0; i < n; i++) scaled[i] = ec[i] * elScale;
        curves.push({ data: scaled, color: `rgba(126,184,255,${alpha.toFixed(2)})`, lineWidth: 0.8 });
      }

      // Sum curve
      const scaledSum = new Float64Array(n);
      for (let i = 0; i < n; i++) scaledSum[i] = sum[i] * sumScale;
      curves.push({
        data: scaledSum, color: '#66ff99', lineWidth: 2.5,
        label: shifted ? 'Sum (reassigned)' : 'Sum (unshifted)'
      });

      yMax = norm ? 1.1 : sumMax * 1.1;
    } else {
      const el = elements[this.selectedElement];
      const elR = Math.sqrt(el.x * el.x + el.y * el.y);
      let elPSF = offAxisPSF2D(grid, el.x, el.y);
      let shiftedPSF = new Float64Array(grid.length);
      const dx = grid[1] - grid[0];
      const shiftBins = (-el.x / 2) / dx;
      for (let i = 0; i < grid.length; i++) {
        const srcIdx = i - shiftBins;
        const lo = Math.floor(srcIdx);
        const frac = srcIdx - lo;
        if (lo >= 0 && lo + 1 < grid.length) {
          shiftedPSF[i] = elPSF[lo] * (1 - frac) + elPSF[lo + 1] * frac;
        }
      }

      if (!norm) {
        // Match the "All" view scaling: compute the same elScale boost
        // so this element appears at the same height as in "All" mode.
        let elAllMax = 0;
        for (const e of elements) {
          const psf = offAxisPSF2D(grid, e.x, e.y);
          for (const v of psf) { if (v > elAllMax) elAllMax = v; }
        }
        // In "All" mode, elements are scaled so tallest reaches 0.35
        // of normalized sum (peak=1). So elScale = 0.35 / elAllMax.
        // Apply that same scale to our curves.
        const elScale = elAllMax > 0 ? 0.35 / elAllMax : 1;
        const scaled = (arr) => {
          const out = new Float64Array(n);
          for (let i = 0; i < n; i++) out[i] = arr[i] * elScale;
          return out;
        };
        elPSF = scaled(elPSF);
        shiftedPSF = scaled(shiftedPSF);
        yMax = 1.1; // same as normalized "All" view
      }

      // Detector element as top-hat at its radial position
      const elSize = 0.35; // element width in AU (matches hex spacing)
      const detTophat = new Float64Array(grid.length);
      // Scale top-hat to ~15% of plot height
      const tophatH = norm ? 1 : yMax * 0.15;
      for (let i = 0; i < grid.length; i++) {
        detTophat[i] = Math.abs(grid[i] - elR) < elSize / 2 ? tophatH : 0;
      }

      curves.push(
        { data: maybeNorm(wf), color: 'rgba(255,255,255,0.2)', lineWidth: 1, dash: [6, 4], label: 'Widefield' },
        { data: detTophat, color: '#ff776644', lineWidth: 1.5, label: `Det. element (r=${elR.toFixed(2)})` },
      );
      if (this._shifted) {
        curves.push(
          { data: maybeNorm(shiftedPSF), color: '#66ff99', lineWidth: 2.5,
            label: `Reassigned (−x/2)` },
        );
      } else {
        curves.push(
          { data: maybeNorm(elPSF), color: '#66ff99', lineWidth: 2.5,
            label: 'Effective PSF' },
        );
      }
    }

    const title = this.selectedElement === -1
      ? (this._shifted ? 'Pixel Reassignment (all elements)' : 'Naive Sum (all elements)')
      : `Detector Element ${this.selectedElement}`;

    drawPlot(ctx, rect, curves, {
      grid, xRange: [-2.5, 2.5], yRange: [0, yMax],
      xLabel: 'Position (AU)', yLabel: norm ? 'Intensity (norm.)' : 'Intensity',
      title,
    });

    // FWHM info
    if (this.selectedElement === -1) {
      // Recompute sum for FWHM (use the already-computed scaledSum approach)
      const n2 = grid.length;
      const dx2 = grid[1] - grid[0];
      const sumForFWHM = new Float64Array(n2);
      for (const el of elements) {
        const psf = offAxisPSF2D(grid, el.x, el.y);
        if (this._shifted) {
          const shiftBins = (-el.x / 2) / dx2;
          for (let i = 0; i < n2; i++) {
            const srcIdx = i - shiftBins;
            const lo = Math.floor(srcIdx);
            const frac = srcIdx - lo;
            if (lo >= 0 && lo + 1 < n2) {
              sumForFWHM[i] += psf[lo] * (1 - frac) + psf[lo + 1] * frac;
            }
          }
        } else {
          for (let i = 0; i < n2; i++) sumForFWHM[i] += psf[i];
        }
      }
      const wfFWHM = measureFWHM(wf, grid);
      const sumFWHM = measureFWHM(sumForFWHM, grid);
      ctx.fillStyle = 'rgba(255,255,255,0.5)';
      ctx.font = '11px system-ui';
      ctx.textAlign = 'left';
      const mode = this._shifted ? 'Reassigned' : 'Unshifted';
      ctx.fillText(
        `Widefield FWHM: ${wfFWHM.toFixed(3)} AU | ${mode} sum FWHM: ${sumFWHM.toFixed(3)} AU | Ratio: ${(wfFWHM / sumFWHM).toFixed(2)}×`,
        rect.x, rect.y + rect.h + 46
      );
    }
  }

  // ── Step 4: Scanning beam build-up ─────────────────────────────

  _renderStep4(ctx, W, H) {
    const grid = this._grid;
    const n = grid.length;
    const dx = grid[1] - grid[0];
    const beamPos = this.beamPos;
    const phSize = this.pinholeSize4;
    const phOffset = this.pinholeOffset4;
    const pinCenter = beamPos + phOffset;
    const norm = this.normalized;

    // Emission PSF scaling from Stokes shift. Airy-disc first zero
    // scales linearly with wavelength, so emission AU = (λ_em/λ_exc) × exc AU.
    // We baseline λ_exc at 488 nm for the display.
    const lambdaExc = 488;
    const emScale = (lambdaExc + this.stokesShift4) / lambdaExc;
    const airyEm = (r) => airyDisk(r / emScale);

    // Excitation beam: Airy disc centered at scan position
    const excBeam = new Float64Array(n);
    for (let i = 0; i < n; i++) excBeam[i] = airyDisk(grid[i] - beamPos);

    // Excitation intensity at the point emitter (x = 0)
    const excAtEmitter = airyDisk(beamPos);

    // Emission PSF: centered at the emitter, scaled by excitation there
    const emissionPSF = new Float64Array(n);
    for (let i = 0; i < n; i++) emissionPSF[i] = airyEm(grid[i]) * excAtEmitter;

    // Static effective PSF with optional pinhole offset D:
    //   I_eff(x_s) = h_exc(x_s) · (h_det ⊛ P)(x_s + D)
    // Derivation: pinhole is centered at x_s + D, so collected signal is
    //   ∫ h_det(x) · tophat(x − (x_s + D), ρ) dx.
    // Since tophat is symmetric this equals (h_det ∗ P)(x_s + D) with P the
    // centered top-hat; we reuse the on-axis convDet and sample it at x_s + D.
    //
    // We compute the *shape* (unscaled by 2·phSize for sub-grid pinholes) so
    // that the normalized-mode display stays meaningful at phSize → 0, where
    // the shape limits to h_exc(x_s) · h_det(x_s + D). The sub-grid scale
    // factor is re-applied below only for unnormalized display, so that the
    // absolute effective PSF still vanishes as phSize → 0.
    const excStatic = airyDiskArray(grid);
    const detStatic = new Float64Array(n);
    for (let i = 0; i < n; i++) detStatic[i] = airyEm(grid[i]);
    // Total emission integral over the simulation grid. We divide the
    // effective PSF by this so the y-axis has a consistent "intensity"
    // interpretation across all curves: effPSF represents excitation ×
    // fraction-of-emission-captured, matching the peak=1 convention used
    // for h_exc and h_det. As phSize → grid extent, effPSF → h_exc.
    const totalDetInt = integrate(detStatic, grid);
    const subGrid = phSize < 2 * dx;
    let convDetShape;
    if (subGrid) {
      convDetShape = detStatic; // limit shape: h_det
    } else {
      const pin = new Float64Array(n);
      for (let i = 0; i < n; i++) pin[i] = tophat(grid[i], phSize);
      convDetShape = convolve1D(detStatic, pin, dx);
    }
    const effPSF = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const srcIdx = (grid[i] + phOffset - grid[0]) / dx;
      const lo = Math.floor(srcIdx);
      const frac = srcIdx - lo;
      let convVal = 0;
      if (lo >= 0 && lo + 1 < n) {
        convVal = convDetShape[lo] * (1 - frac) + convDetShape[lo + 1] * frac;
      } else if (lo === n - 1) {
        convVal = convDetShape[n - 1];
      }
      effPSF[i] = excStatic[i] * convVal / totalDetInt;
    }

    let plotEff = effPSF;
    const yMax = 1.1;
    if (norm) {
      let m = 0;
      for (const v of effPSF) if (v > m) m = v;
      if (m > 0) {
        plotEff = new Float64Array(n);
        for (let i = 0; i < n; i++) plotEff[i] = effPSF[i] / m;
      }
    } else if (subGrid) {
      // Re-apply the 2·phSize factor so unnormalized amplitude collapses to 0
      plotEff = new Float64Array(n);
      for (let i = 0; i < n; i++) plotEff[i] = effPSF[i] * 2 * phSize;
    }

    const rect = { x: 50, y: 30, w: W - 100, h: H - 80 };
    const xMin = -3, xMax = 3;
    const toX = (v) => rect.x + (v - xMin) / (xMax - xMin) * rect.w;
    const toY = (v) => rect.y + rect.h - v / yMax * rect.h;

    // ── Filled region: emission PSF clipped by pinhole window ──
    const pinLeft = Math.max(pinCenter - phSize, xMin);
    const pinRight = Math.min(pinCenter + phSize, xMax);
    if (pinRight > pinLeft && phSize > 0) {
      const interp = (x) => {
        const idx = (x - grid[0]) / dx;
        const lo = Math.floor(idx);
        const frac = idx - lo;
        if (lo < 0) return emissionPSF[0];
        if (lo + 1 >= n) return emissionPSF[n - 1];
        return emissionPSF[lo] * (1 - frac) + emissionPSF[lo + 1] * frac;
      };

      ctx.fillStyle = 'rgba(255,119,102,0.28)';
      ctx.beginPath();
      ctx.moveTo(toX(pinLeft), toY(0));
      ctx.lineTo(toX(pinLeft), toY(interp(pinLeft)));
      for (let i = 0; i < n; i++) {
        if (grid[i] > pinLeft && grid[i] < pinRight) {
          ctx.lineTo(toX(grid[i]), toY(emissionPSF[i]));
        }
      }
      ctx.lineTo(toX(pinRight), toY(interp(pinRight)));
      ctx.lineTo(toX(pinRight), toY(0));
      ctx.closePath();
      ctx.fill();
    }

    // ── Pinhole boundary verticals (faint dashed) ──
    if (phSize > 0) {
      ctx.strokeStyle = 'rgba(255,200,150,0.55)';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      if (pinLeft > xMin) {
        ctx.moveTo(toX(pinLeft), rect.y);
        ctx.lineTo(toX(pinLeft), rect.y + rect.h);
      }
      if (pinRight < xMax) {
        ctx.moveTo(toX(pinRight), rect.y);
        ctx.lineTo(toX(pinRight), rect.y + rect.h);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    const curves = [
      { data: excBeam, color: '#5599ff', lineWidth: 1.5, label: 'Excitation beam' },
      { data: emissionPSF, color: '#ff7766', lineWidth: 2, label: 'Emission PSF' },
      { data: plotEff, color: '#66ff99', lineWidth: 2.5, dash: [6, 4], label: 'Effective PSF' },
      { legendOnly: true, color: 'rgba(255,200,150,0.85)', lineWidth: 1,
        dash: [4, 4], label: 'Pinhole' },
    ];

    drawPlot(ctx, rect, curves, {
      grid, xRange: [xMin, xMax], yRange: [0, yMax],
      xLabel: 'Position (AU)', yLabel: norm ? 'Intensity (norm.)' : 'Intensity',
      title: 'Scanning Confocal: Building the Effective PSF',
    });

    // ── Point emitter marker at x = 0 ──
    ctx.fillStyle = '#ff7766';
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(toX(0), toY(0), 5, 0, 2 * PI);
    ctx.fill();
    ctx.stroke();

    // ── Beam-position marker on effective PSF curve ──
    const beamIdx = Math.round((beamPos - grid[0]) / dx);
    if (beamIdx >= 0 && beamIdx < n) {
      ctx.fillStyle = '#66ff99';
      ctx.strokeStyle = '#fff';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(toX(beamPos), toY(plotEff[beamIdx]), 5, 0, 2 * PI);
      ctx.fill();
      ctx.stroke();
    }
  }

  // ── Reassignment animation ─────────────────────────────────────

  _startReassignmentAnimation() {
    this._reassignProgress = 0;
    const startTime = performance.now();
    const duration = 1500;

    const animate = (now) => {
      this._reassignProgress = Math.min(1, (now - startTime) / duration);
      this._renderFrame();
      if (this._reassignProgress < 1) {
        this._animFrameId = requestAnimationFrame(animate);
      }
    };
    if (this._animFrameId) cancelAnimationFrame(this._animFrameId);
    this._animFrameId = requestAnimationFrame(animate);
  }

  // ── Canvas click handling ──────────────────────────────────────

  _handleCanvasClick(e) {
    if (this.step !== 3 || !this._hexHitData) return;
    const rect = this._canvas.getBoundingClientRect();
    const mx = (e.clientX - rect.left) * (this._canvas.width / rect.width) / (window.devicePixelRatio || 1);
    const my = (e.clientY - rect.top) * (this._canvas.height / rect.height) / (window.devicePixelRatio || 1);

    const { cx, cy, scale, radius, elements, allBtnY } = this._hexHitData;

    // Check "All" button
    if (mx >= cx - 30 && mx <= cx + 30 && my >= allBtnY && my <= allBtnY + 24) {
      this.selectedElement = -1;
      return;
    }

    // Check hex elements
    for (let i = 0; i < elements.length; i++) {
      const el = elements[i];
      const ex = cx + el.x * scale;
      const ey = cy + el.y * scale;
      const dist = Math.sqrt((mx - ex) ** 2 + (my - ey) ** 2);
      if (dist < radius) {
        this.selectedElement = i;
        this._lastSingleElement = i;
        return;
      }
    }
  }

  // ── Event handlers ─────────────────────────────────────────────

  _setStep(step) {
    this.step = step;
    if (step === 2 || step === 3 || step === 4) this.normalized = false;
  }

  _onPinholeSize(e) {
    this.pinholeSize = parseFloat(e.target.value);
  }

  _onPinholeOffset(e) {
    this.pinholeOffset = parseFloat(e.target.value);
  }

  _toggleNormalized() {
    this.normalized = !this.normalized;
  }

  // ── Template ───────────────────────────────────────────────────

  render() {
    const stepDescriptions = {
      1: html`<strong>The pinhole trade-off:</strong> A confocal microscope's PSF is the
        product of its excitation and detection PSFs. Closing the pinhole narrows the
        detection PSF, improving resolution up to <strong>√2×</strong> — but signal
        drops dramatically. Move the slider to see this trade-off.`,
      2: html`<strong>The off-axis surprise:</strong> What happens if the pinhole is
        <em>displaced</em> from the optical axis? The image shifts by exactly
        <strong>half the displacement</strong>, but the resolution is preserved!
        This is the key insight behind ISM.`,
      3: html`<strong>Pixel reassignment:</strong> An array detector (like Airyscan's
        ${this._detectorElements.length} elements) acts as many small off-axis pinholes simultaneously.
        Since each element's image is shifted by a known amount (d/2), we can shift
        them all back and sum — recovering <strong>√2× resolution with full signal</strong>.
        Click detector elements to explore.`,
      4: html`<strong>How the PSF emerges:</strong> A single emitter sits at x=0.
        As the scanning beam moves, it excites the emitter proportionally to its
        intensity there, and the pinhole (moving with the beam) collects only the
        emission that falls inside its window. The <strong>effective PSF</strong>
        is the product of excitation × (emission ∩ pinhole) — try scanning the beam
        and shrinking the pinhole.`,
    };

    return html`
      <div class="container sl-theme-dark">
        <div class="step-nav">
          <button class="step-btn ${this.step === 1 ? 'active' : ''}"
            @click=${() => this._setStep(1)}>1. Pinhole Trade-off</button>
          <button class="step-btn ${this.step === 2 ? 'active' : ''}"
            @click=${() => this._setStep(2)}>2. Off-Axis Effect</button>
          <button class="step-btn ${this.step === 3 ? 'active' : ''}"
            @click=${() => this._setStep(3)}>3. ISM Reconstruction</button>
          <button class="step-btn ${this.step === 4 ? 'active' : ''}"
            @click=${() => this._setStep(4)}>4. PSF Formation</button>
        </div>

        <div class="step-description">${stepDescriptions[this.step]}</div>

        <canvas width="860" height="460"
          @click=${this._handleCanvasClick}></canvas>

        ${this.step === 1 ? this._renderStep1Controls() : ''}
        ${this.step === 2 ? this._renderStep2Controls() : ''}
        ${this.step === 3 ? this._renderStep3Controls() : ''}
        ${this.step === 4 ? this._renderStep4Controls() : ''}
      </div>
    `;
  }

  _renderStep1Controls() {
    return html`
      <div class="controls">
        <div class="control-group">
          <label>Pinhole diameter</label>
          <sl-range min="0.02" max="3" step="0.02"
            value=${this.pinholeSize}
            @sl-input=${this._onPinholeSize}></sl-range>
          <span class="value">${this.pinholeSize.toFixed(2)} AU</span>
        </div>
        <div class="btn-group">
          <button class="${this.normalized ? 'active' : ''}"
            @click=${() => { this.normalized = true; }}>Normalized</button>
          <button class="${!this.normalized ? 'active' : ''}"
            @click=${() => { this.normalized = false; }}>Unnormalized</button>
        </div>
        <div class="btn-group">
          <button class="${this._showConvDet ? 'active' : ''}"
            @click=${() => { this._showConvDet = !this._showConvDet; }}>Show Detection PSF</button>
        </div>
      </div>
      <div class="info-bar">
        <span>Pinhole: <span class="highlight">${this.pinholeSize.toFixed(2)} AU</span></span>
        <span>Relative FWHM: <span class="highlight">${this._getCurrentFWHM()}</span></span>
        <span>Relative signal: <span class="highlight">${this._getCurrentSignal()}</span></span>
      </div>
    `;
  }

  _renderStep2Controls() {
    return html`
      <div class="controls" style="flex-direction: column;">
        <div class="control-group" style="width: 100%;">
          <label>Pinhole offset</label>
          <sl-range min="-2" max="2" step="0.02"
            value=${this.pinholeOffset}
            @sl-input=${this._onPinholeOffset}></sl-range>
          <span class="value">${this.pinholeOffset.toFixed(2)} AU</span>
        </div>
        <div style="display: flex; gap: 20px; align-items: center; width: 100%;">
          <div class="control-group" style="flex: 1;">
            <label>Pinhole size</label>
            <sl-range min="0" max="1" step="0.02"
              value=${this.pinholeSize2}
              @sl-input=${(e) => { this.pinholeSize2 = parseFloat(e.target.value); }}></sl-range>
            <span class="value">${this.pinholeSize2 === 0 ? 'point' : this.pinholeSize2.toFixed(2) + ' AU'}</span>
          </div>
          <div class="btn-group">
            <button class="${this.normalized ? 'active' : ''}"
              @click=${() => { this.normalized = true; }}>Normalized</button>
            <button class="${!this.normalized ? 'active' : ''}"
              @click=${() => { this.normalized = false; }}>Unnormalized</button>
          </div>
        </div>
      </div>
      <div class="info-bar">
        <span>Pinhole offset: <span class="highlight">d = ${this.pinholeOffset.toFixed(2)} AU</span></span>
        <span>Reassignment position: <span class="highlight">d/2 = ${(this.pinholeOffset / 2).toFixed(2)} AU</span></span>
        <span>Resolution preserved ✓</span>
      </div>
    `;
  }

  _renderStep3Controls() {
    const el = this.selectedElement >= 0 ? this._detectorElements[this.selectedElement] : null;
    const elR = el ? Math.sqrt(el.x * el.x + el.y * el.y) : 0;
    return html`
      <div class="controls">
        <div class="btn-group">
          <button class="${this.selectedElement === -1 ? 'active' : ''}"
            @click=${() => { this.selectedElement = -1; }}>All</button>
          <button class="${this.selectedElement !== -1 ? 'active' : ''}"
            @click=${() => { if (this.selectedElement === -1) this.selectedElement = this._lastSingleElement; }}>Single</button>
        </div>
        <div class="btn-group">
          <button class="${!this._shifted ? 'active' : ''}"
            @click=${() => { this._shifted = false; }}>Unshifted</button>
          <button class="${this._shifted ? 'active' : ''}"
            @click=${() => { this._shifted = true; }}>Reassigned</button>
        </div>
        <div class="btn-group">
          <button class="${this.normalized ? 'active' : ''}"
            @click=${() => { this.normalized = true; }}>Normalized</button>
          <button class="${!this.normalized ? 'active' : ''}"
            @click=${() => { this.normalized = false; }}>Unnormalized</button>
        </div>
      </div>
      <div class="info-bar">
        ${this.selectedElement === -1
          ? html`<span>Mode: <span class="highlight">${this._shifted ? 'Pixel reassignment (shifted by d/2)' : 'Naive sum (no shift)'}</span></span>`
          : html`<span>Element: <span class="highlight">${this.selectedElement}</span>
              | Offset r: <span class="highlight">${elR.toFixed(2)} AU</span>
              | Shift: <span class="highlight">${(-elR / 2).toFixed(2)} AU</span></span>`
        }
      </div>
    `;
  }

  _renderStep4Controls() {
    return html`
      <div class="controls" style="flex-direction: column;">
        <div class="control-group" style="width: 100%;">
          <label>Beam position</label>
          <sl-range min="-2" max="2" step="0.02"
            value=${this.beamPos}
            @sl-input=${(e) => { this.beamPos = parseFloat(e.target.value); }}></sl-range>
          <span class="value">${this.beamPos.toFixed(2)} AU</span>
        </div>
        <div style="display: flex; gap: 20px; align-items: center; width: 100%;">
          <div class="control-group" style="flex: 1;">
            <label>Pinhole size</label>
            <sl-range min="0" max="3" step="0.02"
              value=${this.pinholeSize4}
              @sl-input=${(e) => { this.pinholeSize4 = parseFloat(e.target.value); }}></sl-range>
            <span class="value">${this.pinholeSize4 === 0 ? 'point' : this.pinholeSize4.toFixed(2) + ' AU'}</span>
          </div>
          <div class="btn-group">
            <button class="${this.normalized ? 'active' : ''}"
              @click=${() => { this.normalized = true; }}>Normalized</button>
            <button class="${!this.normalized ? 'active' : ''}"
              @click=${() => { this.normalized = false; }}>Unnormalized</button>
          </div>
        </div>
        <div class="control-group" style="width: 100%;">
          <label>Pinhole offset</label>
          <sl-range min="-1" max="1" step="0.02"
            value=${this.pinholeOffset4}
            @sl-input=${(e) => { this.pinholeOffset4 = parseFloat(e.target.value); }}></sl-range>
          <span class="value">${this.pinholeOffset4.toFixed(2)} AU</span>
        </div>
        <div class="control-group" style="width: 100%;">
          <label>Stokes shift</label>
          <sl-range min="0" max="200" step="1"
            value=${this.stokesShift4}
            @sl-input=${(e) => { this.stokesShift4 = parseFloat(e.target.value); }}></sl-range>
          <span class="value">${this.stokesShift4.toFixed(0)} nm</span>
        </div>
      </div>
      <div class="info-bar">
        <span>Beam: <span class="highlight">${this.beamPos.toFixed(2)} AU</span></span>
        <span>Pinhole offset: <span class="highlight">${this.pinholeOffset4.toFixed(2)} AU</span></span>
        <span>Pinhole size: <span class="highlight">${this.pinholeSize4 === 0 ? 'point' : this.pinholeSize4.toFixed(2) + ' AU'}</span></span>
        <span>Stokes shift: <span class="highlight">${this.stokesShift4.toFixed(0)} nm (λ_em/λ_exc = ${((488 + this.stokesShift4) / 488).toFixed(3)})</span></span>
        <span>Excitation at emitter: <span class="highlight">${airyDisk(this.beamPos).toFixed(3)}</span></span>
      </div>
    `;
  }

  _getCurrentFWHM() {
    const sweep = this._pinholeSweep;
    const fwhm = this._fwhmCurve;
    if (!sweep) return '—';
    for (let i = 1; i < sweep.length; i++) {
      if (sweep[i] >= this.pinholeSize) {
        const t = (this.pinholeSize - sweep[i - 1]) / (sweep[i] - sweep[i - 1]);
        return (fwhm[i - 1] + t * (fwhm[i] - fwhm[i - 1])).toFixed(3);
      }
    }
    return fwhm[fwhm.length - 1].toFixed(3);
  }

  _getCurrentSignal() {
    const sweep = this._pinholeSweep;
    const sig = this._signalCurve;
    if (!sweep) return '—';
    for (let i = 1; i < sweep.length; i++) {
      if (sweep[i] >= this.pinholeSize) {
        const t = (this.pinholeSize - sweep[i - 1]) / (sweep[i] - sweep[i - 1]);
        return (sig[i - 1] + t * (sig[i] - sig[i - 1])).toFixed(3);
      }
    }
    return sig[sig.length - 1].toFixed(3);
  }

  disconnectedCallback() {
    super.disconnectedCallback();
    if (this._animFrameId) cancelAnimationFrame(this._animFrameId);
  }
}

customElements.define('ism-tutorial', IsmTutorial);
