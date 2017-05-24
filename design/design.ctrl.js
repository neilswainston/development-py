designApp.controller("designCtrl", function($scope) {
	var self = this;
	
	self.models = {
        selected: null,
        templates: [
        	{
        		type: "design",
        		id: 2,
        		features: []
        	},
            {
        		type: "feature",
        		id: 1
            }
        ],
        designs: [
        	{
        		type: "design",
            	id: 1,
            	features: []
        	}
        ]
    };
	
	self.modelAsJson = angular.toJson(self.models, true)

	$scope.$watch(function() {
		return self.models.designs;
	},               
	function(values) {
		self.modelAsJson = angular.toJson(self.models, true)
	}, true);
});