/** @type {import('tailwindcss').Config} */
module.exports = {
  content: ["./src/**/*.{js,ts,jsx,tsx,mdx}"],
  theme: {
    extend: {
      colors: {
        ink: "#0f1419",
        paper: "#eef1f3",
        accent: "#1a5f4a",
        muted: "#5c6b73",
        panel: "#ffffff",
        "panel-inset": "#f7f9fa",
        "panel-border": "#d5dde2",
      },
      fontFamily: {
        display: ["var(--font-plex-sans)", "ui-sans-serif", "system-ui", "sans-serif"],
        mono: ["var(--font-plex-mono)", "ui-monospace", "monospace"],
      },
    },
  },
  plugins: [],
};
