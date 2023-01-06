/**
 * Implementation of FORBID algorithm
 *
 * FORBID: Fast Overlap Removal By Stochastic Gradient Descent for Graph Drawing,
 * https://github.com/LoannGio/FORBID
 */
import { createFunction, type Algorithm } from 'agora-graph';
import {
  layout,
  loadParams,
  main,
  params,
  type Layout,
  type Parameters,
} from './forbid';

/**
 * Executes the FORBID algorithm for this graph
 *
 * @param {Graph} graph the graph to update
 * @param {object} [options] options to pass to the algorith
 * @param {number} options.padding padding to add between nodes
 *
 * @returns {Result} the updated graph
 */
export const forbid = createFunction(function (
  graph,
  options: Partial<Parameters> = {}
) {
  /* https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L496 */

  const g: Layout = layout(graph.nodes, graph.nodes);

  const result = main(g, loadParams(params, options));
  const X = result[0];
  const S = result[1];

  return {
    graph: {
      nodes: graph.nodes.map((v, i) => ({
        x: X[i * 2],
        y: X[i * 2 + 1],
        width: S[i * 2],
        height: S[i * 2 + 1],
        index: v.index,
        label: v.label,
        meta: v.meta,
        up: v.up,
      })),
      edges: graph.edges,
    },
  };
});

export const FORBIDAlgorithm: Algorithm<Parameters> = {
  name: 'FORBID',
  algorithm: forbid,
};

export default FORBIDAlgorithm;
