# Micro-Active

Interactive web components for microscopy education, designed to be embedded in Moodle lessons.

## Structure

Each component lives in its own directory under `components/`:

```
components/
└── <component-name>/
    ├── <component-name>.js    # Web component (for Moodle embedding)
    └── demo.html              # Local demonstration page
```

## Getting Started

1. Install dependencies:
   ```bash
   npm install
   ```

2. Start a local server:
   ```bash
   npm run dev
   ```

3. Open `http://localhost:8000` in your browser to see all demos.

## Components

- **gaussian-profile**: Interactive Gaussian curve with adjustable width

## Usage in Moodle

To embed a component in Moodle:

1. Upload the component's `.js` file to your Moodle files
2. In your lesson, add an HTML block with:
   ```html
   <script type="module" src="/path/to/component.js"></script>
   <component-name></component-name>
   ```

## Development

Each component is a standalone Lit web component that:
- Bundles all its dependencies via CDN imports
- Can be used independently
- Has a local demo page for testing
