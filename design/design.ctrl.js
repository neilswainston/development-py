designApp.controller("designCtrl", function($scope) {
	var self = this;

	self.models = {
			selected: null,
			templates: [
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
					id: 1,
					name: "promoter",
					temp_params: {
						fixed: true
					}
				},
				{
					type: "feature",
					typ: "http://purl.obolibrary.org/obo/SO_0000139",
					id: 2,
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
					id: 3,
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
					id: 4,
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
				],
				designs: [
					{
						type: "design",
						id: 1,
						name: "Design",
						desc: "Design",
						features: []
					}
					]
	};

	self.modelAsJson = angular.toJson(self.models, true)

	$scope.$watch(function() {
		return self.models;
	},               
	function(values) {
		self.modelAsJson = angular.toJson(self.models, true)
	}, true);
});