import { LitElement, html, css } from 'https://cdn.jsdelivr.net/gh/lit/dist@3/core/lit-core.min.js';

/**
 * Interactive diffraction simulation showing plane wave propagation through
 * single or double slits, computed via Huygens-Fresnel wavelet summation.
 *
 * @element diffraction-slit
 */
class DiffractionSlit extends LitElement {
  static properties = {
    mode: { type: String },
    slitWidth: { type: Number, attribute: 'slit-width' },
    slitSeparation: { type: Number, attribute: 'slit-separation' },
    wavelength: { type: Number },
    renderMode: { type: String, attribute: 'render-mode' },
    numSources: { type: Number, attribute: 'num-sources' },
    waveletCount: { type: Number, attribute: 'wavelet-count' },
    speed: { type: Number },
    highQuality: { type: Boolean, attribute: 'high-quality' },
    isPlaying: { type: Boolean, attribute: 'is-playing' },

    showModeSelector: { type: Boolean, attribute: 'show-mode-selector' },
    showQuality: { type: Boolean, attribute: 'show-quality' },
    showSpeed: { type: Boolean, attribute: 'show-speed' },
    showSlitWidth: { type: Boolean, attribute: 'show-slit-width' },
    showSlitSeparation: { type: Boolean, attribute: 'show-slit-separation' },
    showRenderMode: { type: Boolean, attribute: 'show-render-mode' },
    showWavelength: { type: Boolean, attribute: 'show-wavelength' },
    showPlayButton: { type: Boolean, attribute: 'show-play-button' },
  };

  static styles = css`
    :host {
      display: block;
      font-family: system-ui, -apple-system, sans-serif;
      max-width: 860px;
    }
    .container {
      background: #111122;
      border-radius: 12px;
      padding: 16px;
      box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
    }
    canvas {
      display: block;
      width: 100%;
      border-radius: 8px;
      background: #000;
    }
    .controls {
      display: flex;
      flex-direction: column;
      gap: 10px;
      margin-top: 14px;
      padding: 12px;
      background: rgba(255, 255, 255, 0.04);
      border-radius: 8px;
    }
    .control-row {
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: flex-end;
    }
    .control-group {
      flex: 1;
      min-width: 130px;
    }
    .control-group.fixed {
      flex: 0 0 33%;
      min-width: auto;
    }
    .control-group.narrow {
      flex: 0 0 auto;
      min-width: auto;
    }
    .control-group label {
      display: block;
      color: #9aa;
      font-size: 11px;
      margin-bottom: 4px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }
    .button-group {
      display: flex;
      gap: 3px;
    }
    .button-group button {
      flex: 1;
      padding: 6px 10px;
      border: 1px solid #334;
      background: #1e1e30;
      color: #aab;
      border-radius: 5px;
      cursor: pointer;
      font-size: 12px;
      transition: all 0.15s;
    }
    .button-group button:hover {
      background: #2a2a44;
    }
    .button-group button.active {
      background: #4466dd;
      border-color: #4466dd;
      color: #fff;
    }
    .play-btn {
      padding: 7px 18px;
      border: 1px solid #334;
      background: #1e1e30;
      color: #aab;
      border-radius: 5px;
      cursor: pointer;
      font-size: 13px;
      transition: all 0.15s;
      width: 100%;
    }
    .play-btn:hover { background: #2a2a44; }
    .play-btn.playing { background: #c44; border-color: #c44; color: #fff; }
    .toggle {
      position: relative;
      width: 36px;
      height: 20px;
      background: #333;
      border-radius: 10px;
      cursor: pointer;
      transition: background 0.2s;
      border: none;
      padding: 0;
    }
    .toggle.on { background: #4466dd; }
    .toggle::after {
      content: '';
      position: absolute;
      top: 2px;
      left: 2px;
      width: 16px;
      height: 16px;
      background: #ccc;
      border-radius: 50%;
      transition: transform 0.2s;
    }
    .toggle.on::after { transform: translateX(16px); background: #fff; }
    sl-range {
      --track-color-active: #4466dd;
      --track-color-inactive: #222238;
      --thumb-size: 14px;
      --tooltip-offset: -9999px;
    }
    sl-range::part(tooltip) {
      display: none;
    }
  `;

  constructor() {
    super();
    this.mode = 'double';
    this.slitWidth = 1;
    this.slitSeparation = 5;
    this.wavelength = 20;
    this.renderMode = 'intensity';
    this.numSources = 50;
    this.waveletCount = 7;
    this.speed = 25;
    this.highQuality = true;
    this.isPlaying = true;

    this.showModeSelector = true;
    this.showQuality = true;
    this.showSpeed = true;
    this.showSlitWidth = true;
    this.showSlitSeparation = true;
    this.showRenderMode = true;
    this.showWavelength = true;
    this.showPlayButton = true;

    // Internal state
    this._time = 0;
    this._animationId = null;
    this._recomputeTimer = null;

    // Canvas layout
    this._W = 800;
    this._H = 384;
    this._barrierY = 300;
    this._barrierH = 4;
    this._profileH = 55;

    // Precomputed field arrays (half-resolution grid)
    this._fieldReal = null;
    this._fieldImag = null;
    this._fieldAmp = null;
    this._normScale = 1;
    this._sources = [];
    this._hW = 0;  // reduced-res width
    this._hBY = 0; // half-res barrier height
  }

  async firstUpdated() {
    const imports = [];
    if (!document.querySelector('link[href*="shoelace"]')) {
      const link = document.createElement('link');
      link.rel = 'stylesheet';
      link.href = 'https://cdn.jsdelivr.net/npm/@shoelace-style/shoelace@2.12.0/cdn/themes/dark.css';
      document.head.appendChild(link);
    }
    if (!customElements.get('sl-range'))
      imports.push(import('https://cdn.jsdelivr.net/npm/@shoelace-style/shoelace@2.12.0/cdn/components/range/range.js'));
    if (imports.length) await Promise.all(imports);
    this._canvas = this.renderRoot.querySelector('#wave-canvas');
    this._ctx = this._canvas.getContext('2d');
    this._computeField();
    if (this.isPlaying) this._startAnimation();
    else this._renderFrame();
  }

  disconnectedCallback() {
    this._stopAnimation();
    if (this._recomputeTimer) clearTimeout(this._recomputeTimer);
    super.disconnectedCallback();
  }

  updated(changed) {
    // Prevent slit overlap: separation must be >= slit width
    if (this.mode === 'double' && this.slitSeparation < this.slitWidth) {
      if (changed.has('slitWidth')) this.slitWidth = this.slitSeparation;
      else this.slitSeparation = this.slitWidth;
    }
    const needsRecompute = ['slitWidth', 'slitSeparation', 'mode', 'wavelength', 'numSources', 'highQuality'];
    if (needsRecompute.some(p => changed.has(p)) && this._canvas) {
      this._scheduleRecompute();
    }
    if (changed.has('renderMode') || changed.has('waveletCount')) {
      if (!this.isPlaying && this._canvas) this._renderFrame();
    }
    if (changed.has('isPlaying')) {
      if (this.isPlaying) this._startAnimation();
      else this._stopAnimation();
    }
  }

  _scheduleRecompute() {
    if (this._recomputeTimer) clearTimeout(this._recomputeTimer);
    this._recomputeTimer = setTimeout(() => {
      this._computeField();
      if (!this.isPlaying) this._renderFrame();
    }, 60);
  }

  // --- Slit geometry ---

  _getSlitRanges() {
    const lambda = this.wavelength;
    const cx = this._W / 2;
    const sw = this.slitWidth * lambda;
    if (this.mode === 'single') {
      return [[cx - sw / 2, cx + sw / 2]];
    }
    const sep = this.slitSeparation * lambda;
    return [
      [cx - sep / 2 - sw / 2, cx - sep / 2 + sw / 2],
      [cx + sep / 2 - sw / 2, cx + sep / 2 + sw / 2],
    ];
  }

  // --- Field computation (Huygens-Fresnel wavelet summation) ---

  _computeField() {
    const W = this._W;
    const barrierY = this._barrierY;
    const lambda = this.wavelength;
    const k = (2 * Math.PI) / lambda;
    const res = this.highQuality ? 1 : 2;

    // Adaptive source count: ~12 per wavelength of slit width, min 10
    const N = Math.max(10, Math.round(this.slitWidth * 12));

    // Build source points (canvas x-coords within slits)
    const sources = [];
    const slitRanges = this._getSlitRanges();
    for (const [left, right] of slitRanges) {
      const w = right - left;
      for (let j = 0; j < N; j++) {
        sources.push(left + (j + 0.5) * (w / N));
      }
    }

    const dxSrc = (this.slitWidth * lambda) / N;

    // Compute at reduced resolution
    const hW = Math.ceil(W / res);
    const hBY = Math.ceil(barrierY / res);
    const totalPx = hBY * hW;
    const fR = new Float32Array(totalPx);
    const fI = new Float32Array(totalPx);
    const fA = new Float32Array(totalPx);
    let maxAmp = 0;

    for (let hy = 0; hy < hBY; hy++) {
      const cy = hy * res;
      const dy = barrierY - cy;
      if (dy < 2) continue;
      const dy2 = dy * dy;
      const rowOff = hy * hW;

      for (let hx = 0; hx < hW; hx++) {
        const cx = hx * res;
        let uR = 0, uI = 0;
        for (let s = 0; s < sources.length; s++) {
          const dx = cx - sources[s];
          const r = Math.sqrt(dx * dx + dy2);
          const amp = dy / (r * Math.sqrt(r));
          const phase = k * r;
          uR += amp * Math.cos(phase);
          uI += amp * Math.sin(phase);
        }
        uR *= dxSrc;
        uI *= dxSrc;

        const idx = rowOff + hx;
        fR[idx] = uR;
        fI[idx] = uI;
        const a = Math.sqrt(uR * uR + uI * uI);
        fA[idx] = a;
        if (a > maxAmp) maxAmp = a;
      }
    }

    this._res = res;
    this._hW = hW;
    this._hBY = hBY;
    this._fieldReal = fR;
    this._fieldImag = fI;
    this._fieldAmp = fA;
    this._normScale = maxAmp > 0 ? 1 / maxAmp : 1;
    this._sources = sources;
  }

  // --- Animation ---

  _startAnimation() {
    if (this._animationId) return;
    const tick = (ts) => {
      this._time = ts * 0.001 * this.speed;
      this._renderFrame();
      this._animationId = requestAnimationFrame(tick);
    };
    this._animationId = requestAnimationFrame(tick);
  }

  _stopAnimation() {
    if (this._animationId) {
      cancelAnimationFrame(this._animationId);
      this._animationId = null;
    }
  }

  // --- Rendering ---

  _renderFrame() {
    const ctx = this._ctx;
    if (!ctx || !this._fieldReal) return;

    const W = this._W;
    const H = this._H;
    const lambda = this.wavelength;
    const k = (2 * Math.PI) / lambda;
    const omega = k; // phase velocity = 1 px per time unit
    const t = this._time;
    const barrierY = this._barrierY;
    const barrierH = this._barrierH;
    const barrierBot = barrierY + barrierH;
    const norm = this._normScale;

    const res = this._res;
    const hW = this._hW;
    const cosWt = Math.cos(omega * t);
    const sinWt = Math.sin(omega * t);

    // Build slit membership mask
    const slitRanges = this._getSlitRanges();
    const slitMask = new Uint8Array(W);
    for (const [l, r] of slitRanges) {
      const lo = Math.max(0, Math.floor(l));
      const hi = Math.min(W, Math.ceil(r));
      for (let x = lo; x < hi; x++) slitMask[x] = 1;
    }

    const img = ctx.createImageData(W, H);
    const d = img.data;
    const rm = this.renderMode;

    for (let cy = 0; cy < H; cy++) {
      for (let cx = 0; cx < W; cx++) {
        const pi = (cy * W + cx) * 4;
        let gray = 0;

        if (cy < barrierY) {
          // Diffracted region
          const dy = barrierY - cy;
          if (dy < 2) {
            // Transition row near barrier
            gray = slitMask[cx] ? 80 : 40;
          } else {
            const hx = (cx / res) | 0;
            const hy = (cy / res) | 0;
            const fi = hy * hW + hx;
            const field =
              (this._fieldReal[fi] * cosWt + this._fieldImag[fi] * sinWt) * norm;
            const amp = this._fieldAmp[fi] * norm;

            if (rm === 'intensity') {
              gray = 128 + 127 * clamp(field, -1, 1);
            } else if (rm === 'wavefronts') {
              gray = wavefrontGray(field, amp);
            } else {
              // wavelets: dimmed wave background
              gray = 50 + 30 * clamp(field, -1, 1);
            }
          }
        } else if (cy < barrierBot) {
          // Barrier
          if (slitMask[cx]) {
            // Slit opening
            const f = Math.cos(k * (barrierBot - cy) + omega * t);
            if (rm === 'intensity') gray = 128 + 127 * clamp(f, -1, 1);
            else if (rm === 'wavefronts') gray = wavefrontGray(f, 1);
            else gray = 50 + 30 * clamp(f, -1, 1);
          } else {
            d[pi] = 55; d[pi + 1] = 55; d[pi + 2] = 65; d[pi + 3] = 255;
            continue;
          }
        } else {
          // Plane wave below barrier
          const f = Math.cos(k * (cy - barrierBot) + omega * t);
          if (rm === 'intensity') gray = 128 + 127 * f;
          else if (rm === 'wavefronts') gray = wavefrontGray(f, 1);
          else gray = 50 + 30 * f;
        }

        d[pi] = gray;
        d[pi + 1] = gray;
        d[pi + 2] = gray;
        d[pi + 3] = 255;
      }
    }

    ctx.putImageData(img, 0, 0);

    // Overlay layers
    if (rm === 'wavelets') this._drawWavelets(ctx, k, omega, t);
    this._drawBarrierOverlay(ctx, slitRanges);
    this._drawIntensityProfile(ctx);
  }

  // --- Wavelet arcs overlay ---

  _drawWavelets(ctx, k, omega, t) {
    const lambda = this.wavelength;
    const barrierY = this._barrierY;
    const W = this._W;
    const sources = this._sources;
    if (!sources.length) return;

    // Select evenly-spaced subset
    const perSlit = this.waveletCount;
    const slitRanges = this._getSlitRanges();
    const displaySources = [];
    for (const [left, right] of slitRanges) {
      // find sources within this slit
      const inSlit = sources.filter(sx => sx >= left && sx <= right);
      const n = Math.min(perSlit, inSlit.length);
      if (n <= 0) continue;
      const step = inSlit.length / n;
      for (let i = 0; i < n; i++) {
        displaySources.push(inSlit[Math.floor(i * step + step / 2)]);
      }
    }

    // Phase offset (fractional wavelength)
    const phaseOff = (((omega * t) / (2 * Math.PI)) % 1 + 1) % 1;

    ctx.save();
    // Clip to diffracted region
    ctx.beginPath();
    ctx.rect(0, 0, W, barrierY);
    ctx.clip();

    ctx.strokeStyle = 'rgba(90, 170, 255, 0.45)';
    ctx.lineWidth = 1.2;

    const maxR = Math.sqrt(W * W + barrierY * barrierY);
    for (const sx of displaySources) {
      for (let n = 0; n < maxR / lambda; n++) {
        const r = (n + phaseOff) * lambda;
        if (r < 2) continue;
        ctx.beginPath();
        ctx.arc(sx, barrierY, r, -Math.PI, 0);
        ctx.stroke();
      }
    }
    ctx.restore();
  }

  // --- Barrier overlay ---

  _drawBarrierOverlay(ctx, slitRanges) {
    const W = this._W;
    const by = this._barrierY;
    const bh = this._barrierH;

    ctx.fillStyle = '#383848';
    // Left of first slit
    let prev = 0;
    for (const [l, r] of slitRanges) {
      if (l > prev) ctx.fillRect(prev, by, l - prev, bh);
      prev = r;
    }
    if (prev < W) ctx.fillRect(prev, by, W - prev, bh);

    // Slit edge markers
    ctx.strokeStyle = '#667';
    ctx.lineWidth = 1;
    for (const [l, r] of slitRanges) {
      ctx.strokeRect(l, by, r - l, bh);
    }
  }

  // --- Intensity profile at screen ---

  _drawIntensityProfile(ctx) {
    const W = this._W;
    const ph = this._profileH;
    const norm = this._normScale;
    if (!this._fieldAmp) return;

    // Intensity at the screen (top of diffracted region, hy = 0)
    const hW = this._hW;
    const res = this._res;
    const intensities = new Float32Array(W);
    let maxI = 0;
    for (let cx = 0; cx < W; cx++) {
      const hx = (cx / res) | 0;
      const a = this._fieldAmp[hx]; // row 0 of half-res grid
      const I = a * a * norm * norm;
      intensities[cx] = I;
      if (I > maxI) maxI = I;
    }
    if (maxI === 0) return;

    // Semi-transparent backdrop
    ctx.fillStyle = 'rgba(0, 0, 0, 0.55)';
    ctx.fillRect(0, 0, W, ph);

    // Filled curve
    ctx.beginPath();
    ctx.moveTo(0, ph);
    for (let cx = 0; cx < W; cx++) {
      ctx.lineTo(cx, ph - (intensities[cx] / maxI) * (ph - 6));
    }
    ctx.lineTo(W - 1, ph);
    ctx.closePath();
    ctx.fillStyle = 'rgba(90, 170, 255, 0.25)';
    ctx.fill();
    ctx.strokeStyle = 'rgba(90, 170, 255, 0.7)';
    ctx.lineWidth = 1.5;
    ctx.stroke();

    ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.font = '10px system-ui';
    ctx.fillText('Screen Intensity', 6, 12);
  }

  // --- Template ---

  render() {
    return html`
      <div class="container">
        <canvas id="wave-canvas" width="${this._W}" height="${this._H}"></canvas>
        <div class="controls">
          <div class="control-row">
            ${this.showModeSelector ? html`
              <div class="control-group">
                <label>Slit Configuration</label>
                <div class="button-group">
                  <button class="${this.mode === 'single' ? 'active' : ''}"
                          @click=${() => { this.mode = 'single'; }}>Single</button>
                  <button class="${this.mode === 'double' ? 'active' : ''}"
                          @click=${() => { this.mode = 'double'; }}>Double</button>
                </div>
              </div>
            ` : ''}

            ${this.showSlitWidth ? html`
              <div class="control-group">
                <label>Slit Width: ${this.slitWidth.toFixed(1)}\u03BB</label>
                <sl-range min="0.5" max="10" step="0.1"
                          .value=${this.slitWidth}
                          @sl-input=${(e) => { this.slitWidth = parseFloat(e.target.value); }}>
                </sl-range>
              </div>
            ` : ''}

            ${this.showSlitSeparation && this.mode === 'double' ? html`
              <div class="control-group">
                <label>Separation: ${this.slitSeparation.toFixed(1)}\u03BB</label>
                <sl-range min="0.5" max="25" step="0.1"
                          .value=${this.slitSeparation}
                          @sl-input=${(e) => { this.slitSeparation = parseFloat(e.target.value); }}>
                </sl-range>
              </div>
            ` : ''}

            ${this.showWavelength ? html`
              <div class="control-group">
                <label>Scale (\u03BB = ${this.wavelength}px)</label>
                <sl-range min="8" max="40" step="1"
                          .value=${this.wavelength}
                          @sl-input=${(e) => { this.wavelength = parseInt(e.target.value); }}>
                </sl-range>
              </div>
            ` : ''}
          </div>

          <div class="control-row">
            ${this.showRenderMode ? html`
              <div class="control-group fixed">
                <label>Display Mode</label>
                <div class="button-group">
                  <button class="${this.renderMode === 'intensity' ? 'active' : ''}"
                          @click=${() => { this.renderMode = 'intensity'; }}>Wave</button>
                  <button class="${this.renderMode === 'wavefronts' ? 'active' : ''}"
                          @click=${() => { this.renderMode = 'wavefronts'; }}>Wavefronts</button>
                  <button class="${this.renderMode === 'wavelets' ? 'active' : ''}"
                          @click=${() => { this.renderMode = 'wavelets'; }}>Wavelets</button>
                </div>
              </div>
            ` : ''}

            ${this.renderMode === 'wavelets' ? html`
              <div class="control-group">
                <label>Wavelets: ${this.waveletCount}</label>
                <sl-range min="2" max="20" step="1"
                          .value=${this.waveletCount}
                          @sl-input=${(e) => { this.waveletCount = parseInt(e.target.value); }}>
                </sl-range>
              </div>
            ` : ''}

            ${this.showSpeed ? html`
              <div class="control-group">
                <label>Speed: ${(this.speed / 5).toFixed(0)}\u00D7</label>
                <sl-range min="5" max="100" step="5"
                          .value=${this.speed}
                          @sl-input=${(e) => { this.speed = parseFloat(e.target.value); }}>
                </sl-range>
              </div>
            ` : ''}

            ${this.showQuality ? html`
              <div class="control-group narrow">
                <label>HQ</label>
                <button class="toggle ${this.highQuality ? 'on' : ''}"
                        @click=${() => { this.highQuality = !this.highQuality; }}></button>
              </div>
            ` : ''}

            ${this.showPlayButton ? html`
              <div class="control-group narrow">
                <label>&nbsp;</label>
                <button class="play-btn ${this.isPlaying ? 'playing' : ''}"
                        @click=${() => { this.isPlaying = !this.isPlaying; }}>
                  ${this.isPlaying ? 'Pause' : 'Play'}
                </button>
              </div>
            ` : ''}
          </div>

        </div>
      </div>
    `;
  }
}

// --- Helpers ---

function clamp(v, lo, hi) {
  return v < lo ? lo : v > hi ? hi : v;
}

function wavefrontGray(field, amp) {
  if (amp < 0.005) return 0;
  // normalizedField = cos(phase - wt), ranges -1..1
  const nf = field / amp;
  // Thin bright band at crest
  const crest = Math.max(0, (nf - 0.6) / 0.4);
  // Fade in low-amplitude regions
  const af = Math.min(1, amp * 3);
  return Math.round(crest * af * 255);
}

customElements.define('diffraction-slit', DiffractionSlit);
