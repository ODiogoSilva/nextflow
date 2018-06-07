/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.k8s.client
import javax.net.ssl.HttpsURLConnection

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sClientTest extends Specification {

    def 'should create a request' () {

        given:
        final TOKEN = '8d09d0ds'
        final client = Spy(K8sClient)

        def HTTPS_CONN = Mock(HttpsURLConnection)
        def HTTP_CONN = Mock(HttpURLConnection)

        when:
        client.config.server = 'host.com:443'
        client.config.token = TOKEN
        def resp = client.makeRequest('GET', '/foo/bar')
        then:
        1 * client.createConnection0("https://host.com:443/foo/bar") >> HTTPS_CONN
        1 * client.setupHttpsConn(HTTPS_CONN) >> null
        1 * HTTPS_CONN.setRequestMethod('GET') >> null
        1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer $TOKEN")
        1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json")
        1 * HTTPS_CONN.getResponseCode() >> 200
        1 * HTTPS_CONN.getInputStream() >> { new ByteArrayInputStream('{"field_x":"hello"}'.bytes) }
        resp instanceof K8sResponseApi
        resp.text == '{"field_x":"hello"}'

        when:
        client.config.server = 'http://my-server.com'
        client.config.token = TOKEN
        client.makeRequest('POST', '/foo/bar')
        then:
        1 * client.createConnection0("http://my-server.com/foo/bar") >> HTTP_CONN
        0 * client.setupHttpsConn(_) >> null
        1 * HTTP_CONN.setRequestMethod('POST') >> null
        1 * HTTP_CONN.getResponseCode() >> 401
        1 * HTTP_CONN.getErrorStream() >> { new ByteArrayInputStream('{"field_x":"oops.."}'.bytes) }
        def e = thrown(K8sResponseException)
        e.response.field_x == 'oops..'

    }

    def 'should make a get request' () {

        given:
        def client = Spy(K8sClient)
        when:
        client.get('/foo/bar')
        then:
        1 * client.makeRequest('GET', '/foo/bar') >> null
    }

    def 'should make a post request' () {

        given:
        def client = Spy(K8sClient)
        when:
        client.post('/foo/bar', '{ the: body }')
        then:
        1 * client.makeRequest('POST', '/foo/bar', '{ the: body }') >> null
    }

    def 'should make a delete request' () {

        given:
        def client = Spy(K8sClient)
        when:
        client.delete('/foo/bar', '{ the: body }')
        then:
        1 * client.makeRequest('DELETE', '/foo/bar', '{ the: body }') >> null
    }

    def 'should delete a pod' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"field":"OK"}'

        when:
        result = client.podDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/default/pods/foo",null) >> RESP
        result.field == "OK"

        when:
        client.config.namespace = 'bar'
        result = client.podDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/bar/pods/foo",null) >> RESP
        result.field == "OK"

    }

    def 'should list pods' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.podList(true)
        then:
        1 * client.get('/api/v1/pods') >> RESP
        result.response == 'hello'

        when:
        result = client.podList()
        then:
        1 * client.get('/api/v1/namespaces/default/pods') >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'foo'
        result = client.podList()
        then:
        1 * client.get('/api/v1/namespaces/foo/pods') >> RESP
        result.response == 'hello'
    }

    def 'should list secrets' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.secretesList()
        then:
        1 * client.get("/api/v1/namespaces/default/secrets") >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'pippo'
        result = client.podList()
        then:
        1 * client.get('/api/v1/namespaces/pippo/pods') >> RESP
        result.response == 'hello'

    }

    def 'should describe secret' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.secretDescribe('my-secret')
        then:
        1 * client.get("/api/v1/namespaces/default/secrets/my-secret") >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'paperino'
        result = client.secretDescribe('pluto')
        then:
        1 * client.get('/api/v1/namespaces/paperino/secrets/pluto') >> RESP
        result.response == 'hello'
    }

    def 'should log pod' () {

        given:
        InputStream result
        def client = Spy(K8sClient)
        def STREAM = Mock(InputStream)
        def RESP = Mock(K8sResponseApi)
        RESP.getStream() >> STREAM

        when:
        result = client.podLog('pod-123')
        then:
        1 * client.get('/api/v1/namespaces/default/pods/pod-123/log') >> RESP
        result == STREAM

        when:
        result = client.podLog('pod-123', follow: true)
        then:
        1 * client.get('/api/v1/namespaces/default/pods/pod-123/log?follow=true') >> RESP
        result == STREAM

        when:
        result = client.podLog('pod-123', follow: true, foo:1, bar: 'x')
        then:
        1 * client.get('/api/v1/namespaces/default/pods/pod-123/log?follow=true&foo=1&bar=x') >> RESP
        result == STREAM

    }

    def 'should create config' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'
        def CONFIG = [foo: 'hello', bar:'world']
        def JSON1 = '{"apiVersion":"v1","kind":"ConfigMap","metadata":{"name":"foo","namespace":"default"},"data":{"foo":"hello","bar":"world"}}'
        def JSON2 = '{"apiVersion":"v1","kind":"ConfigMap","metadata":{"name":"foo","namespace":"bar"},"data":{"foo":"hello","bar":"world"}}'

        when:
        result = client.configCreate('foo' , CONFIG)
        then:
        1 * client.post("/api/v1/namespaces/default/configmaps", JSON1) >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'bar'
        result = client.configCreate('foo' , CONFIG)
        then:
        1 * client.post("/api/v1/namespaces/bar/configmaps", JSON2) >> RESP
        result.response == 'done'

    }

    def 'should delete config' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'

        when:
        result = client.configDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/default/configmaps/foo") >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'ns-1'
        result = client.configDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/ns-1/configmaps/foo") >> RESP
        result.response == 'done'

    }

    def 'should delete all configs' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'

        when:
        result = client.configDeleteAll()
        then:
        1 * client.delete("/api/v1/namespaces/default/configmaps") >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'ns-1'
        result = client.configDeleteAll()
        then:
        1 * client.delete("/api/v1/namespaces/ns-1/configmaps") >> RESP
        result.response == 'done'

    }

    def 'should get a pod state' () {

        given:
        def JSON = '''
          {
             "kind": "Pod",
             "apiVersion": "v1",
             "metadata": {
                 "name": "pod-xyz",
                 "namespace": "default",
                 "selfLink": "/api/v1/namespaces/default/pods/pod-xyz/status",
                 "uid": "33390e2c-f84b-11e7-a89d-025000000001",
                 "resourceVersion": "119932",
                 "creationTimestamp": "2018-01-13T10:19:12Z",
                 "labels": {
                     "app": "nextflow"
                 }
             },
             
             
             "status": {
                 "phase": "Succeeded",
                 "conditions": [
                     {
                         "type": "Initialized",
                         "status": "True",
                         "lastProbeTime": null,
                         "lastTransitionTime": "2018-01-13T10:19:12Z",
                         "reason": "PodCompleted"
                     },
                     {
                         "type": "Ready",
                         "status": "False",
                         "lastProbeTime": null,
                         "lastTransitionTime": "2018-01-13T10:19:37Z",
                         "reason": "PodCompleted"
                     },
                     {
                         "type": "PodScheduled",
                         "status": "True",
                         "lastProbeTime": null,
                         "lastTransitionTime": "2018-01-13T10:19:12Z"
                     }
                 ],
                 "hostIP": "192.168.65.3",
                 "podIP": "10.1.0.25",
                 "startTime": "2018-01-13T10:19:12Z",
                 "containerStatuses": [
                     {
                         "name": "pod-xyz",
                         "state": {
                             "terminated": {
                                 "exitCode": 0,
                                 "reason": "Completed",
                                 "startedAt": "2018-01-13T10:19:16Z",
                                 "finishedAt": "2018-01-13T10:19:36Z",
                                 "containerID": "docker://90d447d4a2518642b12c8979474aa2bbb5fe8e96ed9e5caf3c979bb3b751519a"
                             }
                         },
                         "lastState": {

                         },
                         "ready": false,
                         "restartCount": 0,
                         "image": "debian:latest",
                         "imageID": "docker-pullable://debian@sha256:0a5fcee6f52d5170f557ee2447d7a10a5bdcf715dd7f0250be0b678c556a501b",
                         "containerID": "docker://90d447d4a2518642b12c8979474aa2bbb5fe8e96ed9e5caf3c979bb3b751519a"
                     }
                 ],
                 "qosClass": "BestEffort"
             }
         }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'pod-xyz'

        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)

        result == [terminated: [exitCode:0,
                                reason: 'Completed',
                                startedAt: "2018-01-13T10:19:16Z",
                                finishedAt: "2018-01-13T10:19:36Z",
                                containerID: "docker://90d447d4a2518642b12c8979474aa2bbb5fe8e96ed9e5caf3c979bb3b751519a"]]


    }

    def 'should return undetermined status' () {
        given:
        def JSON = '''
             {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-eb853c8010b8e173b23d8d15489d1a31",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-eb853c8010b8e173b23d8d15489d1a31/status",
                     "uid": "5ccdeb23-4f69-11e8-89b1-fa163e31bb09",
                     "resourceVersion": "2847182",
                     "creationTimestamp": "2018-05-04T07:04:18Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "markDuplicates",
                         "runName": "grave-jones",
                         "sessionId": "uuid-f51cd941-c21c-447b-86ca-eaebafa5ad9b",
                         "taskName": "markDuplicates_22028_2_118_1AlignedByCoord.out"
                     }
                 },
              
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-04T07:04:18Z"
                         }
                     ],
                     "qosClass": "Guaranteed"
                 }
             }        
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-eb853c8010b8e173b23d8d15489d1a31'
        
        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)

        result == [:]
    }

    def 'should fail to get pod state' () {

        given:
        def client = Spy(K8sClient)
        final POD_NAME = 'pod-xyz'
        final STATE = [foo:1, bar:2]

        when:
        def e
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([:])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s invalid pod status (missing container status)')

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: []]])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s invalid pod status (missing container status)')

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: [ [name: 'foo'] ]]])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s invalid pod status (name does not match)')

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: [ [name: POD_NAME] ]]])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s invalid pod status (missing state object)')

        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: [ [name: POD_NAME, state: STATE] ]]])
        result == STATE
    }

    def 'client should throw an exception when container status returns ErrImagePull' () {
        given:
        def JSON = '''
             {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-c6e49ee9ebc79486f774a47924a743d7/status",
                     "uid": "18954b48-5811-11e8-8e71-025000000001",
                     "resourceVersion": "3615842",
                     "creationTimestamp": "2018-05-15T07:25:08Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "sayHello",
                         "runName": "lethal-jones",
                         "sessionId": "uuid-b6243a3d-5c07-44f4-97eb-9de499747800",
                         "taskName": "sayHello_2"
                     }
                 },
                 "spec": {
 
                 },
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "Initialized",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         },
                         {
                             "type": "Ready",
                             "status": "False",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z",
                             "reason": "ContainersNotReady",
                             "message": "containers with unready status: [nf-c6e49ee9ebc79486f774a47924a743d7]"
                         },
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         }
                     ],
                     "hostIP": "192.168.65.3",
                     "podIP": "10.1.3.173",
                     "startTime": "2018-05-15T07:25:08Z",
                     "containerStatuses": [
                         {
                             "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                             "state": {
                                 "waiting": {
                                     "reason": "ErrImagePull",
                                     "message": "rpc error: code = Unknown desc = Error response from daemon: pull access denied for nextflow/foo, repository does not exist or may require 'docker login'"
                                 }
                             },
                             "lastState": {
                                 
                             },
                             "ready": false,
                             "restartCount": 0,
                             "image": "nextflow/foo",
                             "imageID": ""
                         }
                     ],
                     "qosClass": "BestEffort"
                 }
             }
'''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-c6e49ee9ebc79486f774a47924a743d7'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(PodUnschedulableException)
        e.message == "K8s pod image cannot be pulled -- rpc error: code = Unknown desc = Error response from daemon: pull access denied for nextflow/foo, repository does not exist or may require 'docker login'"
    }

    def 'client should throw an exception when container status returns ImagePullBackOff' () {
        def JSON = '''
            {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-c6e49ee9ebc79486f774a47924a743d7/status",
                     "uid": "18954b48-5811-11e8-8e71-025000000001",
                     "resourceVersion": "3615866",
                     "creationTimestamp": "2018-05-15T07:25:08Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "sayHello",
                         "runName": "lethal-jones",
                         "sessionId": "uuid-b6243a3d-5c07-44f4-97eb-9de499747800",
                         "taskName": "sayHello_2"
                     }
                 },
                 "spec": { },
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "Initialized",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         },
                         {
                             "type": "Ready",
                             "status": "False",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z",
                             "reason": "ContainersNotReady",
                             "message": "containers with unready status: [nf-c6e49ee9ebc79486f774a47924a743d7]"
                         },
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         }
                     ],
                     "hostIP": "192.168.65.3",
                     "podIP": "10.1.3.173",
                     "startTime": "2018-05-15T07:25:08Z",
                     "containerStatuses": [
                         {
                             "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                             "state": {
                                 "waiting": {
                                     "reason": "ImagePullBackOff",
                                     "message": "Back-off pulling image \\"nextflow/foo\\""
                                 }
                             },
                             "lastState": {
                                 
                             },
                             "ready": false,
                             "restartCount": 0,
                             "image": "nextflow/foo",
                             "imageID": ""
                         }
                     ],
                     "qosClass": "BestEffort"
                 }
             }
'''
        def client = Spy(K8sClient)
        final POD_NAME = 'nf-c6e49ee9ebc79486f774a47924a743d7'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(PodUnschedulableException)
        e.message == "K8s pod image cannot be pulled -- Back-off pulling image \"nextflow/foo\""

    }

    def 'client should throw an exception when pod is Unschedulable' () {

        given:
        def JSON = '''
             {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-52b131a30381b0b88926a9c12e5b1ff1",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-52b131a30381b0b88926a9c12e5b1ff1/status",
                     "uid": "7f17fea0-4e47-11e8-89b1-fa163e31bb09",
                     "resourceVersion": "2674905",
                     "creationTimestamp": "2018-05-02T20:29:21Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "star",
                         "runName": "pedantic-legentil",
                         "sessionId": "uuid-bb3f1a1f-ad9a-4e5f-a2f1-46a4d3fbd2a1",
                         "taskName": "star_22028_2_118_1"
                     }
                 },
                 "spec": { },
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "PodScheduled",
                             "status": "False",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-02T20:29:21Z",
                             "reason": "Unschedulable",
                             "message": "0/4 nodes are available: 4 Insufficient cpu, 4 Insufficient memory."
                         }
                     ],
                     "qosClass": "Guaranteed"
                 }
             }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-52b131a30381b0b88926a9c12e5b1ff1'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(PodUnschedulableException)
        e.message == 'K8s pod cannot be scheduled -- 0/4 nodes are available: 4 Insufficient cpu, 4 Insufficient memory.'
    }
}