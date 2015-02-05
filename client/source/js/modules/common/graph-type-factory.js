define(['angular', 'underscore'], function (angular, _) {
    'use strict';

    return angular.module('app.common.graph-type',[])
    .factory('graphTypeFactory', ['CONFIG', function (CONFIG) {

      /**
       * Iterate through all the types & disable the once where data is missing.
       */
      var enableAnnualCostOptions = function (types, graphData) {
        return _(types.financialAnnualCosts).each(function(entry) {
          if (graphData && graphData.costann &&
            graphData.costann.existing &&
            graphData.costann.existing[entry.id] &&
            graphData.costann.existing[entry.id].legend &&
            graphData.costann.existing[entry.id].legend.length) {
              entry.disabled = false;
            }
          });
      };

      var resetAnnualCostOptions = function(types) {
        return _(types.financialAnnualCosts).each(function(entry) {
          entry.disabled = true;
        });
      };

      return {
          types: angular.copy(CONFIG.GRAPH_TYPES),
          enableAnnualCostOptions: enableAnnualCostOptions,
          resetAnnualCostOptions: resetAnnualCostOptions
      };
    }]);
});
