define(['angular'], function (module) {
  'use strict';

  return angular.module('app.more-than', [])
    /**
     * Adds a more than validator to the parsers of an input element.
     *
     * The input element is only valid if the actual value is larger
     * than the provided value through the directive.
     */
    .directive('moreThan', [function () {
      return {
        require: 'ngModel',
        scope: {
          moreThan: '=',
        },
        link: function(scope, element, attrs, ctrl) {

          var updateValidity = function (highValue, lowValue) {
            var isValid = highValue > lowValue;
            ctrl.$setValidity('moreThan', isValid);
          };

          var validator = function (value) {
            updateValidity(value, scope.moreThan);
            return value;
          };

          // updated validity as soon as the input value changes
          ctrl.$parsers.unshift(validator);
          ctrl.$formatters.unshift(validator);

          // updated validity also when the value to compare changes
          scope.$watch('moreThan', function(moreThanValue) {
            updateValidity(ctrl.$modelValue, moreThanValue);
          });
        }
      };
    }]);
});
