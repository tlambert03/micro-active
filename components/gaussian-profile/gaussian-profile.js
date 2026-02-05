import { LitElement, html, css } from 'https://cdn.jsdelivr.net/gh/lit/dist@3/core/lit-core.min.js';

/**
 * A web component that displays an interactive Airy disk pattern
 * demonstrating diffraction-limited resolution in microscopy.
 *
 * @element gaussian-profile
 */
class GaussianProfile extends LitElement {
  static properties = {
    na: { type: Number },
    wavelength: { type: Number },
    distance: { type: Number },
    pixelSize: { type: Number },
    pixelOffset: { type: Number },
    noiseLevel: { type: Number },
    isPlaying: { type: Boolean }
  };

  static styles = css`
    :host {
      display: block;
      font-family: system-ui, -apple-system, sans-serif;
      padding: 20px;
      max-width: 1000px;
    }

    .container {
      background: white;
      border-radius: 8px;
      padding: 20px;
      box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
    }

    h2 {
      margin: 0 0 20px 0;
      color: #333;
      font-size: 1.5em;
    }

    .controls {
      margin-bottom: 20px;
      display: flex;
      flex-direction: column;
      gap: 15px;
    }

    .control-group {
      display: flex;
      align-items: center;
      gap: 10px;
    }

    label {
      font-weight: 500;
      color: #555;
      min-width: 120px;
    }

    sl-range {
      flex: 1;
    }

    .value-display {
      min-width: 80px;
      text-align: right;
      font-weight: 600;
      color: #2c3e50;
    }

    .button-row {
      margin-top: 15px;
      display: flex;
      gap: 10px;
      width: 100%;
    }

    .play-button, .criterion-button {
      padding: 10px 20px;
      font-size: 0.95em;
      font-weight: 600;
      border: none;
      border-radius: 6px;
      cursor: pointer;
      color: white;
      transition: background 0.2s;
    }

    .play-button {
      background: #3b82f6;
      flex: 2;
      min-height: 44px;
      display: flex;
      align-items: center;
      justify-content: center;
    }

    .play-button:hover {
      background: #2563eb;
    }

    .play-button.playing {
      background: #ef4444;
    }

    .play-button.playing:hover {
      background: #dc2626;
    }

    .criterion-button {
      font-size: 0.9em;
      padding: 8px 16px;
      flex: 1;
      background: #94a3b8;
      opacity: 0.6;
      transition: all 0.2s;
    }

    .criterion-button:hover {
      background: #64748b;
    }

    .criterion-button.active {
      background: #10b981;
      opacity: 1;
      box-shadow: 0 0 0 3px rgba(16, 185, 129, 0.3);
    }

    .criterion-button.active:hover {
      background: #059669;
    }

    .info-box {
      background: #f0f7ff;
      border-left: 4px solid #3b82f6;
      padding: 12px 15px;
      margin-bottom: 20px;
      border-radius: 4px;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }

    .info-box p {
      margin: 0;
      color: #1e40af;
      font-size: 0.95em;
      line-height: 1.5;
    }

    .sampling-info {
      display: flex;
      align-items: center;
      gap: 10px;
    }

    .nyquist-button {
      padding: 6px 12px;
      font-size: 0.85em;
      font-weight: 600;
      border: none;
      border-radius: 4px;
      cursor: pointer;
      background: #0ea5e9;
      color: white;
      transition: background 0.2s;
    }

    .nyquist-button:hover {
      background: #0284c7;
    }

    .chart-container {
      position: relative;
      height: 450px;
      margin-top: 20px;
    }

    canvas {
      max-height: 100%;
    }
  `;

  constructor() {
    super();
    this.na = 1.0;           // Numerical aperture
    this.wavelength = 550;    // Wavelength in nanometers
    this.distance = 0.5;      // Distance between objects in micrometers
    this.pixelSize = 0;       // Pixel size in micrometers (0 = continuous)
    this.pixelOffset = 0;     // Pixel grid offset in micrometers
    this.noiseLevel = 0;      // RMS noise level
    this.isPlaying = false;   // Animation state
    this.chart = null;
    this.animationInterval = null;
  }

  firstUpdated() {
    this.loadDependencies();
  }

  async loadDependencies() {
    // Load Web Awesome components
    if (!customElements.get('sl-range')) {
      const link = document.createElement('link');
      link.rel = 'stylesheet';
      link.href = 'https://cdn.jsdelivr.net/npm/@shoelace-style/shoelace@2.12.0/cdn/themes/light.css';
      document.head.appendChild(link);

      await import('https://cdn.jsdelivr.net/npm/@shoelace-style/shoelace@2.12.0/cdn/components/range/range.js');
    }

    // Load Chart.js
    if (typeof Chart === 'undefined') {
      await import('https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.js');
    }

    // Initialize chart after dependencies are loaded
    this.initChart();
  }

  /**
   * Bessel function of the first kind, order 1
   * Uses series expansion for accurate calculation
   */
  besselJ1(x) {
    if (x === 0) return 0;

    const ax = Math.abs(x);
    if (ax < 8.0) {
      // Series expansion for small x
      const y = x * x;
      let result = x * (72362614232.0 + y * (-7895059235.0 + y *
        (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y *
        (-30.16036606))))));
      result /= (144725228442.0 + y * (2300535178.0 + y *
        (18583304.74 + y * (99447.43394 + y * (376.9991397 + y)))));
      return result;
    } else {
      // Asymptotic expansion for large x
      const z = 8.0 / ax;
      const y = z * z;
      const xx = ax - 2.356194491;
      const p1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 +
        y * (0.2457520174e-5 + y * (-0.240337019e-6))));
      const q1 = 0.04687499995 + y * (-0.2002690873e-3 +
        y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
      const result = Math.sqrt(0.636619772 / ax) *
        (Math.cos(xx) * p1 - z * Math.sin(xx) * q1);
      return x < 0 ? -result : result;
    }
  }

  /**
   * Calculate Airy disk intensity pattern
   * Formula: I(r) = I₀ [2J₁(x)/x]²
   * where x = (2π/λ) * NA * r
   * center: position of the Airy disk center
   */
  airyDisk(r, wavelength, na, center = 0) {
    const rFromCenter = Math.abs(r - center);

    if (rFromCenter === 0) return 1.0;

    // Convert wavelength from nm to μm
    const lambda = wavelength / 1000;

    // Calculate the argument for the Bessel function
    const x = (2 * Math.PI / lambda) * na * rFromCenter;

    // Airy pattern: [2*J1(x)/x]²
    const j1 = this.besselJ1(x);
    const intensity = Math.pow(2 * j1 / x, 2);

    return intensity;
  }

  /**
   * Generate Airy disk data for two objects separated by distance d
   */
  generateAiryData(wavelength, na, distance) {
    // Fixed range: -2 to +2 micrometers
    const rMin = -2.0;
    const rMax = 2.0;
    const numPoints = 300;

    // Object positions: centered at ±d/2
    const pos1 = -distance / 2;
    const pos2 = distance / 2;

    const data1 = [];
    const data2 = [];
    const dataSum = [];

    for (let i = 0; i < numPoints; i++) {
      const rVal = rMin + (rMax - rMin) * i / (numPoints - 1);

      const i1 = this.airyDisk(rVal, wavelength, na, pos1);
      const i2 = this.airyDisk(rVal, wavelength, na, pos2);
      const iSum = i1 + i2;

      data1.push({ x: rVal, y: i1 });
      data2.push({ x: rVal, y: i2 });
      dataSum.push({ x: rVal, y: iSum });
    }

    return { data1, data2, dataSum };
  }

  /**
   * Generate Gaussian random noise using Box-Muller transform
   */
  gaussianRandom(mean = 0, stdDev = 1) {
    const u1 = Math.random();
    const u2 = Math.random();
    const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    return z0 * stdDev + mean;
  }

  /**
   * Discretize continuous data into pixels
   */
  discretizeData(continuousData, pixelSize, pixelOffset, noiseLevel = 0) {
    if (pixelSize <= 0) return [];

    const discretized = [];
    const rMin = -2.0;
    const rMax = 2.0;

    // Calculate pixel boundaries
    // Start from offset and create pixels of size pixelSize
    let pixelStart = -pixelOffset;

    // Adjust to cover the full range
    while (pixelStart > rMin) {
      pixelStart -= pixelSize;
    }

    // Process each pixel
    while (pixelStart < rMax) {
      const pixelEnd = pixelStart + pixelSize;
      const pixelCenter = (pixelStart + pixelEnd) / 2;

      // Find all continuous points that fall within this pixel
      const pointsInPixel = continuousData.filter(
        p => p.x >= pixelStart && p.x < pixelEnd
      );

      if (pointsInPixel.length > 0) {
        // Calculate average intensity for this pixel
        let avgIntensity = pointsInPixel.reduce((sum, p) => sum + p.y, 0) / pointsInPixel.length;

        // Add Gaussian noise
        if (noiseLevel > 0) {
          avgIntensity += this.gaussianRandom(0, noiseLevel);
        }

        // Create step: add points at both edges with same y value
        discretized.push({ x: pixelStart, y: avgIntensity });
        discretized.push({ x: pixelEnd, y: avgIntensity });
      }

      pixelStart = pixelEnd;
    }

    return discretized;
  }

  initChart() {
    const canvas = this.shadowRoot.querySelector('canvas');
    if (!canvas) return;

    const ctx = canvas.getContext('2d');

    // Generate initial data
    const { data1, data2, dataSum } = this.generateAiryData(
      this.wavelength,
      this.na,
      this.distance
    );

    const discretizedData = this.discretizeData(
      dataSum,
      this.pixelSize,
      this.pixelOffset,
      this.noiseLevel
    );

    this.chart = new Chart(ctx, {
      type: 'line',
      data: {
        datasets: [
          {
            label: 'Object 1',
            data: data1,
            borderColor: 'rgba(59, 130, 246, 0.6)',
            backgroundColor: 'rgba(59, 130, 246, 0.05)',
            borderWidth: 1.5,
            borderDash: [5, 5],
            fill: false,
            tension: 0,
            pointRadius: 0
          },
          {
            label: 'Object 2',
            data: data2,
            borderColor: 'rgba(239, 68, 68, 0.6)',
            backgroundColor: 'rgba(239, 68, 68, 0.05)',
            borderWidth: 1.5,
            borderDash: [5, 5],
            fill: false,
            tension: 0,
            pointRadius: 0
          },
          {
            label: 'Combined Image',
            data: dataSum,
            borderColor: 'rgb(16, 185, 129)',
            backgroundColor: 'rgba(16, 185, 129, 0.1)',
            borderWidth: 2.5,
            fill: true,
            tension: 0,
            pointRadius: 0,
            hidden: this.pixelSize > 0
          },
          {
            label: 'Pixelated Image',
            data: discretizedData,
            borderColor: 'rgb(147, 51, 234)',
            backgroundColor: 'rgba(147, 51, 234, 0.15)',
            borderWidth: 2.5,
            fill: true,
            tension: 0,
            pointRadius: 0,
            stepped: false,
            hidden: this.pixelSize === 0
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: {
            display: true,
            position: 'top',
            labels: {
              usePointStyle: true,
              pointStyle: 'line',
              padding: 15
            }
          },
          title: {
            display: true,
            text: 'Airy Disk - Point Spread Function',
            font: { size: 16 }
          },
          tooltip: {
            mode: 'index',
            intersect: false,
            callbacks: {
              label: (context) => {
                return `Intensity: ${context.parsed.y.toFixed(4)}`;
              },
              title: (context) => {
                return `Position: ${context.parsed.x.toFixed(3)} μm`;
              }
            }
          }
        },
        scales: {
          x: {
            type: 'linear',
            title: {
              display: true,
              text: 'Distance (μm)'
            },
            min: -2.0,
            max: 2.0,
            grid: {
              display: false
            },
            border: {
              display: false
            },
            ticks: {
              maxTicksLimit: 11
            }
          },
          y: {
            display: false,
            min: -0.4,
            max: 2.0,
            grid: {
              display: true,
              color: (context) => {
                // Make the y=0 line solid and dark
                return context.tick.value === 0 ? '#000' : 'transparent';
              },
              lineWidth: (context) => {
                return context.tick.value === 0 ? 2 : 0;
              },
              drawTicks: false
            }
          }
        },
        animation: false
      }
    });
  }

  handleNAChange(e) {
    this.na = parseFloat(e.target.value);
    this.updateChart();
  }

  handleWavelengthChange(e) {
    this.wavelength = parseFloat(e.target.value);
    this.updateChart();
  }

  handleDistanceChange(e) {
    this.distance = parseFloat(e.target.value);
    this.updateChart();
  }

  handlePixelSizeChange(e) {
    this.pixelSize = parseFloat(e.target.value);
    this.updateChart();
  }

  handlePixelOffsetChange(e) {
    this.pixelOffset = parseFloat(e.target.value);
    this.updateChart();
  }

  handleNoiseLevelChange(e) {
    this.noiseLevel = parseFloat(e.target.value);
    this.updateChart();
  }

  togglePlayback() {
    this.isPlaying = !this.isPlaying;

    if (this.isPlaying) {
      // Start animation at 15 fps (66.67ms interval)
      this.animationInterval = setInterval(() => {
        this.updateChart();
      }, 1000 / 15);
    } else {
      // Stop animation
      if (this.animationInterval) {
        clearInterval(this.animationInterval);
        this.animationInterval = null;
      }
    }
  }

  disconnectedCallback() {
    super.disconnectedCallback();
    // Clean up animation interval when component is removed
    if (this.animationInterval) {
      clearInterval(this.animationInterval);
    }
  }

  setRayleighCriterion() {
    // d = 0.61 * λ / NA
    const lambda = this.wavelength / 1000; // Convert nm to μm
    this.distance = 0.61 * lambda / this.na;
    this.updateChart();
  }

  setAbbeCriterion() {
    // d = λ / (2 * NA) = 0.5 * λ / NA
    const lambda = this.wavelength / 1000; // Convert nm to μm
    this.distance = 0.5 * lambda / this.na;
    this.updateChart();
  }

  setSparrowCriterion() {
    // d ≈ 0.47 * λ / NA (where combined profile becomes flat at center)
    const lambda = this.wavelength / 1000; // Convert nm to μm
    this.distance = 0.47 * lambda / this.na;
    this.updateChart();
  }

  isNearCriterion(criterionValue, tolerance = 0.01) {
    return Math.abs(this.distance - criterionValue) < tolerance;
  }

  calculateRayleigh() {
    const lambda = this.wavelength / 1000;
    return 0.61 * lambda / this.na;
  }

  calculateAbbe() {
    const lambda = this.wavelength / 1000;
    return 0.5 * lambda / this.na;
  }

  calculateSparrow() {
    const lambda = this.wavelength / 1000;
    return 0.47 * lambda / this.na;
  }

  calculateAiryRadius() {
    const lambda = this.wavelength / 1000;
    return 0.61 * lambda / this.na;
  }

  getSamplingRate() {
    if (this.pixelSize === 0) return 'Continuous';
    const airyRadius = this.calculateAiryRadius();
    const rate = airyRadius / this.pixelSize;
    return rate.toFixed(2);
  }

  setNyquistSampling() {
    const airyRadius = this.calculateAiryRadius();
    this.pixelSize = airyRadius / 2; // 2 pixels per Airy radius
    this.updateChart();
  }

  updateChart() {
    if (!this.chart) return;

    const { data1, data2, dataSum } = this.generateAiryData(
      this.wavelength,
      this.na,
      this.distance
    );
    const discretizedData = this.discretizeData(
      dataSum,
      this.pixelSize,
      this.pixelOffset,
      this.noiseLevel
    );

    this.chart.data.datasets[0].data = data1;
    this.chart.data.datasets[1].data = data2;
    this.chart.data.datasets[2].data = dataSum;
    this.chart.data.datasets[2].hidden = this.pixelSize > 0;
    this.chart.data.datasets[3].data = discretizedData;
    this.chart.data.datasets[3].hidden = this.pixelSize === 0;
    this.chart.update();
  }

  getRayleighCriterion() {
    // Rayleigh criterion: r = 0.61λ/NA (in μm)
    return (0.61 * this.wavelength / 1000 / this.na).toFixed(3);
  }

  getAbbeDiffraction() {
    // Abbe diffraction limit: d = λ/(2NA) (in μm)
    return (this.wavelength / 1000 / (2 * this.na)).toFixed(3);
  }

  render() {
    return html`
      <div class="container">
        <h2>Airy Disk Pattern</h2>

        <div class="info-box">
          <p>
            <strong>Rayleigh Criterion:</strong> ${this.getRayleighCriterion()} μm |
            <strong>Abbe Limit:</strong> ${this.getAbbeDiffraction()} μm
          </p>
          <div class="sampling-info">
            <p>
              <strong>Sampling Rate:</strong> ${this.getSamplingRate()}${this.pixelSize > 0 ? ' px/Airy' : ''}
            </p>
            <button
              class="nyquist-button"
              @click="${this.setNyquistSampling}"
            >
              Nyquist
            </button>
          </div>
        </div>

        <div class="controls">
          <div class="control-group">
            <label for="na">Numerical Aperture:</label>
            <sl-range
              id="na"
              min="0.1"
              max="1.4"
              step="0.05"
              value="${this.na}"
              @sl-input="${this.handleNAChange}"
            ></sl-range>
            <span class="value-display">${this.na.toFixed(2)}</span>
          </div>

          <div class="control-group">
            <label for="wavelength">Wavelength (nm):</label>
            <sl-range
              id="wavelength"
              min="400"
              max="700"
              step="10"
              value="${this.wavelength}"
              @sl-input="${this.handleWavelengthChange}"
            ></sl-range>
            <span class="value-display">${this.wavelength} nm</span>
          </div>

          <div class="control-group">
            <label for="distance">Separation (μm):</label>
            <sl-range
              id="distance"
              min="0"
              max="2.0"
              step="0.01"
              value="${this.distance}"
              @sl-input="${this.handleDistanceChange}"
            ></sl-range>
            <span class="value-display">${this.distance.toFixed(2)} μm</span>
          </div>

          <div class="control-group">
            <label for="pixelSize">Pixel Size (μm):</label>
            <sl-range
              id="pixelSize"
              min="0"
              max="1.0"
              step="0.01"
              value="${this.pixelSize}"
              @sl-input="${this.handlePixelSizeChange}"
            ></sl-range>
            <span class="value-display">${this.pixelSize.toFixed(2)} μm</span>
          </div>

          <div class="control-group">
            <label for="pixelOffset">Pixel Offset (μm):</label>
            <sl-range
              id="pixelOffset"
              min="0"
              max="${this.pixelSize > 0 ? (this.pixelSize / 2).toFixed(3) : 0.5}"
              step="0.001"
              value="${this.pixelOffset}"
              @sl-input="${this.handlePixelOffsetChange}"
            ></sl-range>
            <span class="value-display">${this.pixelOffset.toFixed(3)} μm</span>
          </div>

          <div class="control-group">
            <label for="noiseLevel">Noise Level (RMS):</label>
            <sl-range
              id="noiseLevel"
              min="0"
              max="0.5"
              step="0.01"
              value="${this.noiseLevel}"
              @sl-input="${this.handleNoiseLevelChange}"
            ></sl-range>
            <span class="value-display">${this.noiseLevel.toFixed(2)}</span>
          </div>

          <div class="button-row">
            <button
              class="play-button ${this.isPlaying ? 'playing' : ''}"
              @click="${this.togglePlayback}"
            >
              ${this.isPlaying ? '⏸ Stop' : '▶ Play'}
            </button>
            <button
              class="criterion-button ${this.isNearCriterion(this.calculateRayleigh()) ? 'active' : ''}"
              @click="${this.setRayleighCriterion}"
            >
              Rayleigh Criterion
            </button>
            <button
              class="criterion-button ${this.isNearCriterion(this.calculateAbbe()) ? 'active' : ''}"
              @click="${this.setAbbeCriterion}"
            >
              Abbe Limit
            </button>
            <button
              class="criterion-button ${this.isNearCriterion(this.calculateSparrow()) ? 'active' : ''}"
              @click="${this.setSparrowCriterion}"
            >
              Sparrow Criterion
            </button>
          </div>
        </div>

        <div class="chart-container">
          <canvas></canvas>
        </div>
      </div>
    `;
  }
}

customElements.define('gaussian-profile', GaussianProfile);
