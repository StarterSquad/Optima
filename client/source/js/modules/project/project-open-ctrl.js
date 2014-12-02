// ProjectOpenController deals with loading and removing projects

define(['./module', 'angular', 'underscore'], function (module, angular, _) {
  'use strict';

  module.controller('ProjectOpenController', function ($scope, $http, activeProject, localStorage, projects, modalService) {

    $scope.projects = _.map(projects.projects, function(project){
      project.creation_time = Date.parse(project.creation_time);
      project.data_upload_time = Date.parse(project.data_upload_time);
      return project;
    });

    /**
     * Opens an existing project using `name`
     *
     * Alerts the user if it cannot do it.
     */
    $scope.open = function (name) {
      $http.get('/api/project/open/' + name)
        .success(function (response) {
          if (response && response.status === 'NOK') {
            alert(response.reason);
            return;
          }
          activeProject.setValue(name);
        });
    };

    /**
     * Regenerates workbook for the given project `name`
     * Alerts the user if it cannot do it.
     *
     */
    $scope.workbook = function (name) {
      // read that this is the universal method which should work everywhere in
      // http://stackoverflow.com/questions/24080018/download-file-from-a-webapi-method-using-angularjs
      window.open('/api/project/workbook/' + name, '_blank', '');  
    };

    /**
     * Removes the project
     *
     * If the removed project is the active one it will reset it alerts the user
     * in case of failure.
     */
    var removeNoQuestionsAsked = function (name, index) {
      $http.delete('/api/project/delete/' + name)
        .success(function (response) {
          if (response && response.status === 'NOK') {
            alert(response.reason);
            return;
          }

          $scope.projects.splice(index, 1);

          if (activeProject.name === name) {
            activeProject.setValue('');
          }
        })
        .error(function () {
          alert('Could not remove the project');
        });
    };

    /**
     * Opens a modal window to ask the user for confirmation to remove the project and
     * removes the project if the user confirms. 
     * Closes it without further action otherwise.
     */
    $scope.remove = function ($event, name, index) {
      if ($event) { $event.preventDefault(); }
      modalService.confirm(
        function (){console.log('onAccepted');}, 
        function (){console.log('onCancel');}, 
        'I can haz the modal', 
        'Oh hai!'
      );


      // theModal.model = {
      //     title: 'Remove project',
      //     message: ,
      //     confirmText: 'Yes, remove this project',
      //     cancelText: 'No',
      //     onAccepted: function (){
      //       console.log('onConfirm');
      //       removeNoQuestionsAsked(name, index)},
      //     onRejected: function (){
      //       console.log('onCancel');
      //       return false}
      //   };
    };
  });

});
