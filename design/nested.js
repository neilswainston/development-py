/**
 * The controller doesn't do much more than setting the initial data model
 */
angular.module("demo").controller("NestedListsDemoController", function($scope) {

    $scope.models = {
        selected: null,
        templates: [
            {type: "item", typeName: "Component", id: 2},
            {type: "container", typeName: "Design", id: 1, columns: [[]]}
        ],
        dropzone: [
            {
                "type": "container",
                "id": 1,
                "columns": [[]]
            },
        ]
    };

    $scope.$watch('models.dropzone', function(model) {
        $scope.modelAsJson = angular.toJson(model, true);
    }, true);

});

