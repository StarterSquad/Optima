define(['./module', 'angular', 'underscore'], function (module, angular, _) {
  'use strict';

  module.controller('ProjectCreateOrEditController', function ($scope, $state, $modal,
    $timeout, activeProject, parametersResponse, defaultsResponse, info,
    UserManager, modalService,projects) {

    $scope.allProjectNames = _(projects.projects).map(function(project){return project.name});

    $scope.projectExists = function(){
      var exists = isEditMode()? false:_($scope.allProjectNames).contains($scope.projectParams.name);
      $scope.CreateOrEditProjectForm.ProjectName.$setValidity("projectExists", !exists);
      return exists;
    };

    $scope.projectParams = {
      name: ''
    };
    $scope.editParams = {
      isEdit: false
    };
    $scope.projectInfo = info;

    var availableParameters = parametersResponse.data.parameters;
    var availableDefaults = defaultsResponse.data;

    $scope.submit = "Create project & Optima template";
    $scope.populations = availableDefaults.populations;
    $scope.programs = availableDefaults.programs;
    $scope.categories = availableDefaults.categories;

    function isEditMode(){
      return $state.current.name == "project.edit";
    }

    if (isEditMode()) {
      // change submit button name
      $scope.submit = "Save project & Optima template";

      $scope.editParams.isEdit = true;
      $scope.editParams.canUpdate = true;
      $scope.oldProjectName =  $scope.projectInfo.name;

      if (activeProject.isSet()) {
        $scope.projectParams.name = $scope.oldProjectName;

        $scope.projectParams.datastart = $scope.projectInfo.dataStart;
        $scope.projectParams.dataend = $scope.projectInfo.dataEnd;
        $scope.projectParams.econ_dataend = $scope.projectInfo.projectionEndYear;
      }
      _($scope.populations).each(function(population){
        var source = _.findWhere($scope.projectInfo.populations, { short_name: population.short_name });
        if (source) {
          population.active = true;
          _.extend(population, source);
        }
      });
      _($scope.programs).each(function(program){
        var source = _.findWhere($scope.projectInfo.programs, { short_name: program.short_name });
        if (source) {
          program.active = true;
          _(program).extend(angular.copy(source));
          _(program.parameters).each(function(parameter){
            parameter.active = true;
          });
        }
      });
    }

    // Helper function to open a population modal
    var openPopulationModal = function (population) {
      return $modal.open({
        templateUrl: 'js/modules/project/create-population-modal.html',
        controller: 'ProjectCreatePopulationModalController',
        resolve: {
          population: function () {
            return population;
          }
        }
      });
    };

    /*
     * Filters programs by category
     */
    $scope.filterPrograms = function(category) {
      return _($scope.programs).filter(function (item) {
          return item.category==category.category;
        });
    };

    /*
     * Creates a new population and opens a modal for editing.
     *
     * The entry is only pushed to the list of populations if editing in the modal
     * ended with a successful save.
     */
    $scope.openAddPopulationModal = function ($event) {
      if ($event) {
        $event.preventDefault();
      }
      var population = {};

      return openPopulationModal(population).result.then(
        function (newPopulation) {
          $scope.populations.push(newPopulation);
        }
      );
    };

    /*
     * Opens a modal for editing an existing population.
     */
    $scope.openEditPopulationModal = function ($event, population) {
      if ($event) {
        $event.preventDefault();
      }

      return openPopulationModal(population).result.then(
        function (newPopulation) {
          _(population).extend(newPopulation);
        }
      );
    };

    /*
     * Makes a copy of an existing population and opens a modal for editing.
     *
     * The entry is only pushed to the list of populations if editing in the
     * modal ended with a successful save.
     */
    $scope.copyPopulationAndOpenModal = function ($event, existingPopulation) {
      if ($event) {
        $event.preventDefault();
      }
      var population = angular.copy(existingPopulation);

      return openPopulationModal(population).result.then(
        function (newPopulation) {
          $scope.populations.push(newPopulation);
        }
      );
    };

    // Helper function to open a program modal
    var openProgramModal = function (program) {
      return $modal.open({
        templateUrl: 'js/modules/project/create-program-modal.html',
        controller: 'ProjectCreateProgramModalController',
        resolve: {
          program: function () {
            return program;
          },
          availableParameters: function () {
            return availableParameters;
          },
          populations: function () {
            var activePopulations = toCleanArray($scope.populations);
            return activePopulations;
          }
        }
      });
    };

    /*
     * Creates a new program and opens a modal for editing.
     *
     * The entry is only pushed to the list of programs if editing in the modal
     * ended with a successful save.
     */
    $scope.openAddProgramModal = function ($event) {
      if ($event) {
        $event.preventDefault();
      }
      var program = {};

      return openProgramModal(program).result.then(
        function (newProgram) {
          $scope.programs.push(newProgram);
        }
      );
    };

    /*
     * Opens a modal for editing an existing program.
     */
    $scope.openEditProgramModal = function ($event, program) {
      if ($event) {
        $event.preventDefault();
      }

      return openProgramModal(program).result.then(
        function (newProgram) {
          _(program).extend(newProgram);
        }
      );
    };

    /*
     * Makes a copy of an existing program and opens a modal for editing.
     *
     * The entry is only pushed to the list of programs if editing in the modal
     * ended with a successful save.
     */
    $scope.copyProgram = function ($event, existingProgram) {
      if ($event) {
        $event.preventDefault();
      }
      var program = angular.copy(existingProgram);

      return openProgramModal(program).result.then(
        function (newProgram) {
          $scope.programs.push(newProgram);
        }
      );
    };

    /*
    * Returns true of the two provided arrays are identic
    */
    var areEqualArrays = function(arrayOne, arrayTwo) {
      return _(arrayOne).every(function(element, index) {
        return element === arrayTwo[index];
      });
    };

    /*
     * Returns a collection of entries where all non-active antries are filtered
     * out and the active attribute is removed from each of these entries.
     */
    var toCleanArray = function (collection) {
      return _(collection).chain()
        .where({ active: true })
        .map(function (item) {
          var cl = _(item).omit(['active', '$$hashKey']);
          if (cl.parameters) {
            cl.parameters = _(cl.parameters).chain()
              .where({ active: true })
              .map(function (param) {
                return _(param).omit(['active', '$$hashKey']);
              })
              .value();
            if (cl.parameters === 0) delete cl.parameters;
          }
          return cl;
        })
        .value();
    };

    /**
     * Returns a list of error messages for invalid programs.
     *
     * A program is not valid if it contains a parameter with a population not
     * contained in the available populations
     */
    var findInvalidProgramParameters = function(selectedPopulations, selectedPrograms) {
      var avialablePopulationOptions = _(selectedPopulations).pluck('short_name');
      avialablePopulationOptions.push('ALL_POPULATIONS');

      var parameterErrors = [];

      _(selectedPrograms).each(function(program) {
        _(program.parameters).each(function(parameter) {
          if (!_(avialablePopulationOptions).contains(parameter.value.pops[0])) {
            parameterErrors.push("Program \"" + program.name + "\" contains a parameter with a population not selected for this Project.");
          }
        });
      });

      return _(parameterErrors).uniq().map(function(message) {
        return {message: message};
      });
    };

    /*
     * Returns the provide programs with every "ALL_POPULATIONS" entry replaced
     * by the selected populations.
     *
     * Example: ['ALL_POPULATIONS'] -> ["FSW","CSW","MSM","PWID","CHILD","INF"]
     */
    var insertSelectedPopulations = function (programs, selectedPopulations) {
      var shortPopulationNames = _(selectedPopulations).pluck('short_name');
      return _(programs).map(function(program) {
        program.parameters = _(program.parameters).map(function(entry) {
          if (entry.value.pops[0] === "ALL_POPULATIONS") {
            entry.value.pops = shortPopulationNames;
          }
          return entry;
        });
        return program;
      });
    };

    $scope.prepareCreateOrEditForm = function () {
      if ($scope.CreateOrEditProjectForm.$invalid) {
        modalService.informError([{message: 'Please fill in all the required project fields'}]);
        return false;
      }

      var selectedPrograms = toCleanArray($scope.programs);
      var selectedPopulations = toCleanArray($scope.populations);
      var parameterErrors = findInvalidProgramParameters(selectedPopulations, selectedPrograms);

      if (!_(parameterErrors).isEmpty()) {
        modalService.informError(parameterErrors);
        return false;
      }

      if ( $state.current.name == "project.edit" ) {
        if ( !angular.equals( selectedPopulations,$scope.projectInfo.populations ) ||
             !angular.equals( selectedPrograms,$scope.projectInfo.programs ) ) {
          $scope.editParams.canUpdate = $scope.editParams.canUpdate && selectedPopulations.length == $scope.projectInfo.populations.length;
          $scope.editParams.canUpdate = $scope.editParams.canUpdate && selectedPrograms.length == $scope.projectInfo.programs.length;
          var message = 'You have made changes to populations and programs. All existing data will be lost. Would you like to continue?';
          if ($scope.editParams.canUpdate) {
            message = 'You have changed some program or population parameters. Your original data can be reapplied, but you will have to redo the calibration and analysis. Would you like to continue?';
          }
          modalService.confirm(
            function (){ continueSubmitForm( selectedPrograms, selectedPopulations ); },
            function (){},
            'Yes, save this project',
            'No',
            message,
            'Save Project?'
          );
        } else {
          var message = 'No parameters have been changed. Do you intend to reload the original data and start from scratch?';
          modalService.confirm(
            function (){ continueSubmitForm( selectedPrograms, selectedPopulations ); },
            function (){},
            'Yes, reload this project',
            'No',
            message,
            'Reload project?'
          );
        }
      } else {
        continueSubmitForm( selectedPrograms, selectedPopulations );
      }
    };

    // handle another function to continue to submit form
    // since the confirm modal is async and doesn't wait for user's response
    var continueSubmitForm = function( selectedPrograms, selectedPopulations ) {
      var params = _($scope.projectParams).omit('name');
      params.populations = selectedPopulations;
      params.programs = insertSelectedPopulations(selectedPrograms, selectedPopulations);

      $scope.formAction = '/api/project/create/' + $scope.projectParams.name;
      $scope.formParams = JSON.stringify(params);

      // according to documentation it should have been working without this line, but no cigar
      // https://docs.angularjs.org/api/ng/directive/ngSubmit
      document.getElementById('createForm').action = $scope.formAction;
      document.getElementById('params').value = $scope.formParams;
      document.getElementById('edit_params').value = JSON.stringify($scope.editParams);
      document.getElementById('createForm').submit();

      // update active project
      activeProject.setActiveProjectFor($scope.projectParams.name, UserManager.data);

      // Hack to wait for the project to be created.
      // There is not easy way to intercept the completion of the form submission...
      $timeout(function () {
        $state.go('home');
      }, 3000);

      return true;
    }

  });

});
