define(['angular', 'underscore'], function (angular, _) {
  'use strict';

  return angular.module('app.common.error-messages', []).directive('errorMessages', function () {
    return {
      restrict: 'EA',
      replace: false,
      transclude: false,
      scope: {
        for: '@',
        rules: "="
      },
      require: '^form',
      template: '<div ng-if="errorMessages().length>0" class="error-hint"><div ng-repeat="message in errorMessages()">{{message}}</div></div>',
      link: function ($scope, $elem, $attrs, form) {
        var errorMessages = {
          'min': 'The minimum value must be <%= min %>',
          'max': 'The maximum value can be <%= max %>',
          'required': 'The field <%= name %> is required',
          'greaterThan': 'The value must be greater than <%= greaterThan %>.'
        };
        $scope.errorMessages = function () {
          if (form && form[$scope.for].$dirty) {
            return _.compact(_(form[$scope.for].$error).map(function (e, key) {
              if (e) {
                var template = {};
                if ($scope.rules[key]) {
                  /*
                    If the key is 'required', and the rules object contains 'name' property,
                    show the field name in the error message
                    otherwise just say that the field is required
                  */
                  if ($scope.rules.name && key === 'required') {
                    template.name = $scope.rules.name;
                  }
                  else {
                    template[key] = $scope.rules[key];
                  }
                }
                return _.template(errorMessages[key], template);
              }
            }));
          }
        };
      }
    }
  });
});
