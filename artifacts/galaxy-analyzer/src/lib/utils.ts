import { type ClassValue, clsx } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}

export function formatScientific(num: number): string {
  if (num === 0) return "0";
  if (Math.abs(num) < 0.01 || Math.abs(num) > 9999) {
    return num.toExponential(2);
  }
  return num.toFixed(2);
}
