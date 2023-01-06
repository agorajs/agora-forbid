import { type Point, type Box } from 'agora-graph';
type Coord = Point;
type Size = Box;
export type Layout = {
    pos: Coord[];
    sizes: Size[];
    N: number;
};
export declare function layout(pos: Coord[], sizes: Size[]): {
    pos: Point[];
    sizes: Box[];
    N: number;
};
export type Parameters = {
    K: number;
    ALPHA: number;
    MINIMUM_MOVEMENT: number;
    MAX_ITER: number;
    MAX_PASSES: number;
    SCALE_STEP: number;
    distance_fn: (ci: Coord, cj: Coord, si: Size, sj: Size, intersec_width: number, intersec_height: number) => number;
    delta: number;
    eps: number;
    seed: number;
    PRIME: boolean;
};
export declare const params: Parameters;
/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L289 */
export declare function loadParams(defaultParams: Parameters, params: Partial<Parameters>): Parameters;
/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L496 */
export declare function main(g: Layout, options: Parameters): number[][];
export {};
//# sourceMappingURL=forbid.d.ts.map