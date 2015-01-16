define(['./module'], function (module) {
  'use strict';

  module.controller('TypeSelectorController', ['$scope', '$state',
    function ($scope, $state) {

        /**
        * @description Adds a parameter on scope that will indicate if something should be visible depending on the current state
        * @param part name of the parameter that should be set on $scope
        * @param viewNames names of the views in which this parameter should return true
        */
        function showPartInViews(part, viewNames){
            $scope[part+"IsVisible"] = _(viewNames).contains($state.current.name);
        }

        $scope.$on('$stateChangeSuccess',function(){
          showPartInViews('typeSelector', [
            'model.view',
            'analysis.scenarios',
            'analysis.optimization'
          ]);

          showPartInViews('stackedCheckbox',[
            'model.view'
          ]);
        });
    }
  ]);
});
