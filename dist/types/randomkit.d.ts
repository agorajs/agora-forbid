export type rk_state = {
    lcg: () => number;
};
/**
 * Initialize the RNG state using the given seed.
 */
export declare function rk_seed(seed: number): rk_state;
/**
 * Returns a random unsigned long between 0 and ULONG_MAX inclusive
 */
export declare function rk_interval(max: number, state: rk_state): number;
//# sourceMappingURL=randomkit.d.ts.map