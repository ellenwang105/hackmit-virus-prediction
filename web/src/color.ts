import { scaleLinear } from "d3-scale";

/**
 * Epitopes are ~5% of residues, so a linear 0-1 scale would leave nearly
 * everything the same pale colour. The ramp runs over sqrt(score) instead;
 * the legend ticks are still real probabilities.
 */
const ramp = scaleLinear<string>()
  .domain([0, 0.3, 0.6, 1])
  .range(["#ccd5e1", "#f4cf5c", "#f0803c", "#a4161a"])
  .clamp(true);

const LUT_SIZE = 256;
const lutCss: string[] = [];
const lutInt: number[] = [];
for (let k = 0; k < LUT_SIZE; k++) {
  const css = ramp(k / (LUT_SIZE - 1));
  const [r, g, b] = (css.match(/\d+/g) ?? ["0", "0", "0"]).map(Number);
  lutCss.push(css);
  lutInt.push((r << 16) | (g << 8) | b);
}
const bucket = (score: number) =>
  Math.round(Math.sqrt(Math.min(1, Math.max(0, score))) * (LUT_SIZE - 1));

export const scoreCss = (score: number) => lutCss[bucket(score)];
/** 0xRRGGBB, the form 3Dmol colour functions want */
export const scoreInt = (score: number) => lutInt[bucket(score)];

export const LEGEND_TICKS = [0, 0.05, 0.25, 0.5, 1];
export const legendGradient = () =>
  `linear-gradient(to right, ${[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    .map((t) => `${ramp(t)} ${t * 100}%`)
    .join(", ")})`;

export const TRUTH_COLOR = "#0f9d9a";
export const SELECT_COLOR = "#7c3aed";
export const GHOST_COLOR = "#aab3c0";
export const GLYCAN_COLOR = "#b7a3d3";
export const hex = (css: string) => parseInt(css.slice(1), 16);

export const SITE_COLOR: Record<string, string> = {
  A: "#4e79a7",
  B: "#f28e2b",
  C: "#59a14f",
  D: "#b07aa1",
  E: "#9c755f",
};
export const REGION_COLOR = { head: "#7f8ea3", stem: "#c2cad6" } as const;
export const FUNCTION_COLOR = { rbs: "#d1495b", fusion: "#3b7ea1" } as const;
