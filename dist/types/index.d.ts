/**
 * Implementation of FORBID algorithm
 *
 * FORBID: Fast Overlap Removal By Stochastic Gradient Descent for Graph Drawing,
 * https://github.com/LoannGio/FORBID
 */
import { type Algorithm } from 'agora-graph';
import { type Parameters } from './forbid';
/**
 * Executes the FORBID algorithm for this graph
 *
 * @param {Graph} graph the graph to update
 * @param {object} [options] options to pass to the algorith
 * @param {number} options.padding padding to add between nodes
 *
 * @returns {Result} the updated graph
 */
export declare const forbid: import("agora-graph").Function<Partial<Parameters>>;
export declare const FORBIDAlgorithm: Algorithm<Parameters>;
export default FORBIDAlgorithm;
//# sourceMappingURL=index.d.ts.map