import { randomInt, randomLcg } from 'd3-random';

export type rk_state = {
  lcg: () => number;
};

/**
 * Initialize the RNG state using the given seed.
 */
export function rk_seed(seed: number): rk_state {
  return { lcg: randomLcg(seed) };
}

/**
 * Returns a random unsigned long between 0 and ULONG_MAX inclusive
 */
export function rk_interval(max: number, state: rk_state) {
  return randomInt.source(state.lcg)(max + 1)();
}
