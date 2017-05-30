designApp.controller("designCtrl", ["$scope", "DesignService", function($scope, DesignService) {
	var self = this;
	
	self.templates = [
		{
			typ: "http://purl.obolibrary.org/obo/SO_0001416",
			name: "5' flanking region",
			temp_params: {
				fixed: true
			}
		},
		{
			typ: "http://purl.obolibrary.org/obo/SO_0000143",
			name: "assembly component",
			temp_params: {
				fixed: true
			}
		},
		{
			type: "feature",
			typ: "http://purl.obolibrary.org/obo/SO_0000167",
			name: "promoter",
			temp_params: {
				fixed: true
			}
		},
		{
			type: "feature",
			typ: "http://purl.obolibrary.org/obo/SO_0000139",
			name:"ribosome entry site",
			end: 60,
			parameters: {
				"TIR target": 15000
			},
			temp_params: {
				fixed: false
			}
		},
		{
			type: "feature",
			typ: "http://purl.obolibrary.org/obo/SO_0000316",
			name: "coding sequence",
			options: [
				{
					typ: "http://purl.obolibrary.org/obo/SO_0000316",
					name: "coding sequence",
					temp_params: {
						fixed: false
					}
				}
				]
		},
		{
			type: "feature",
			typ: "http://purl.obolibrary.org/obo/SO_0000141",
			name: "terminator",
			temp_params: {
				fixed: true
			}
		},
		{
			typ: "http://purl.obolibrary.org/obo/SO_0000449",
			end: 100,
			name: "random region",
			temp_params: {
				fixed: false
			}
		},
		{
			typ: "http://purl.obolibrary.org/obo/SO_0001417",
			name: "3' flanking region",
			temp_params: {
				fixed: true
			}
		}
	];
	
	self.query = function() {
		return DesignService.query;
	};
	
	self.selected = function() {
		return DesignService.selected;
	};
	
	self.setSelected = function(selected) {
		return DesignService.setSelected(selected);
	};
	
	self.addDesign = function() {
		self.query().designs.push({
			type: "design",
			name: "Design",
			desc: "Design",
			features: []
		});
	};
	
	self.removeDesign = function(index) {
		self.query().designs.splice(index, 1);
	};

	self.queryJson = angular.toJson({selected: self.selected(), query: self.query()}, true);
	
	$scope.$watch(function() {
		return self.selected();
	},               
	function(values) {
		self.queryJson = angular.toJson({selected: self.selected(), query: self.query()}, true)
	}, true);
	
	$scope.$watch(function() {
		return self.query();
	},               
	function(values) {
		self.queryJson = angular.toJson({selected: self.selected(), query: self.query()}, true)
	}, true);
	
	// Initialise UI:
	self.addDesign();
}]);