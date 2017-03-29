define(['./module', 'sha224/sha224'], function (module, SHA224) {
  'use strict';

  return module.controller('RegisterController', function ($scope, $window, userApi) {

    $scope.error = false;

    $scope.register = function () {
      $scope.$broadcast('form-input-check-validity');

      if ($scope.RegisterForm.$invalid) {
        return;
      }

      $scope.error = false;

      var hashed_password = SHA224($scope.password).toString();

      userApi.create({
        username: $scope.username,
        password: hashed_password,
        displayName: $scope.displayName,
        email: $scope.email
      },
        // success
        function (response) {
          if (response.username) {
            // success
            $window.location = './#/login';
          }
        },
        // error
        function (error) {
          $scope.error = error.data.reason;
          switch(error.status){
            case 409: // conflict: will be sent if the email already exists
              // show css error tick to email field
              $scope.RegisterForm.email.$invalid = true;
              $scope.RegisterForm.email.$valid = false;
              $scope.$broadcast('form-input-check-validity');
              break;
            case 400:
              break;
            default:
              $scope.error = 'HTTP status code:' + error.status;
          }
        }
      );
    };

  });

});
