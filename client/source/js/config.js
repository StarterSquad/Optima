﻿/**
 * Defines constants for application
 */
define(['angular'], function (angular) {
  'use strict';
  return angular.module('app.constants', [])
    .constant('CONFIG', {
      GRAPH_MARGINS: {
        top: 20,
        right: 20,
        bottom: 45,
        left: 70
      },
      GRAPH_TYPES: {
        plotUncertainties: false,
        timeVaryingOptimizations: false,
        population: [
          { id: 'prev', name: 'HIV prevalence', byPopulation: true, total: false},
          { id: 'plhiv', name:'Number of PLHIV', byPopulation:false, total:true, stacked: true},
          { id: 'daly', name: 'HIV-related DALYs', byPopulation: false, total: true, stacked: true },
          { id: 'death', name: 'AIDS-related deaths', byPopulation: false, total: true, stacked: true },
          { id: 'inci', name: 'New HIV infections', byPopulation: false, total: true, stacked: true },
          { id: 'dx', name: 'New HIV diagnoses', byPopulation: false, total: true, stacked: true },
          { id: 'tx1', name: 'People on 1st-line treatment', byPopulation: false, total: true, stacked: true },
          { id: 'tx2', name: 'People on 2nd-line treatment', byPopulation: false, total: true, stacked: true }
        ],
        financial: [
          { id: 'total', name: 'Total HIV-related financial commitments', annual: true, cumulative: true  },
          { id: 'existing', name: 'Financial commitments for existing PLHIV', annual: true, cumulative: true },
          { id: 'future', name: 'Financial commitments for future PLHIV', annual: true, cumulative: true }
        ],
        financialAnnualCosts: [
          {id:'total', name:'Total amount', disabled: false},
          {id:'gdp', name:'Proportion of GDP', disabled: false},
          {id:'revenue', name:'Proportion of government revenue', disabled: false},
          {id:'govtexpend', name:'Proportion of government expenditure', disabled: false},
          {id:'totalhealth', name:'Proportion of total health expenditure', disabled: false},
          {id:'domestichealth', name:'Proportion of domestic health expenditure', disabled:false}
        ],
        annualCost: 'total'
      }
    });
});
