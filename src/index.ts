import * as THREE from 'three';
import * as dat from 'dat.gui';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';
import Stats from 'three/examples/jsm/libs/stats.module.js';

interface ModelData {
  points: number[][];
  salinity: number[];
  temperature: number[];
  velocity: number[][];
}
interface Data {
  salinity: number[];
  temperature: number[];
  velocity: number[];
}
type ColorValue = [number, number, number];
const hi: ColorValue = [180 / 255, 6 / 255, 38 / 255];
const mi: ColorValue = [218 / 255, 220 / 255, 223 / 255];
const lo: ColorValue = [60 / 255, 78 / 255, 194 / 255];
const colorMap = {
  hi, mi, lo
}
let currentDimension: string = 'salinity';
//scene
let sphereScene: THREE.Scene;
let planeScene: THREE.Scene;
let currentScene: THREE.Scene;

let geometry = new THREE.BufferGeometry();
let renderer: THREE.WebGLRenderer;
let scene: THREE.Scene;
let camera: THREE.PerspectiveCamera;
let light: THREE.AmbientLight;
let modelData: ModelData;
let stats: Stats;
let controls: OrbitControls;
let currentPoints: THREE.Points;
let data: Data = {
  salinity: [],
  temperature: [],
  velocity: []
}
let minSalinity: number;
let maxSalinity: number;
let minTemperature: number;
let maxTemperature: number;
let minVelocity: number;
let maxVelocity: number;

document.body.onload = function () {
  //init GUI
  initGui();
  initRender();
  initScene();
  initCamera();
  initLight();
  initModel();
  initControls();
  initStats();
  animate();

  // Create scenes
  const sphereScene = new THREE.Scene();
  const planeScene = new THREE.Scene();

  // Set the initial scene


  // Create sphere and plane geometries, materials, and meshes
  const sphereGeometry = new THREE.SphereGeometry(50, 32, 32);
  const planeGeometry = new THREE.PlaneGeometry(100, 100);
  const sphereMaterial = new THREE.PointsMaterial({ size: 0.3, vertexColors: true });
  const planeMaterial = new THREE.PointsMaterial({ size: 0.3, vertexColors: true });
  const sphereMesh = new THREE.Points(sphereGeometry, sphereMaterial);
  const planeMesh = new THREE.Points(planeGeometry, planeMaterial);
  // Add meshes to scenes
  sphereScene.add(sphereMesh);
  planeScene.add(planeMesh);

  // create scene change button
  const sceneChangeButton = document.createElement('button');
  sceneChangeButton.innerText = 'Change Scene';
  sceneChangeButton.style.position = 'absolute';
  sceneChangeButton.style.top = '10px';
  sceneChangeButton.style.left = '10px';
  
  sceneChangeButton.addEventListener('click', () => {
    if(currentScene === sphereScene) {
      onSceneChange('plane');
      // currentScene = planeScene;
      // scene.remove(sphereScene);
      // scene.add(planeScene);
    } else {
      onSceneChange('sphere');
      // currentScene = sphereScene;
      // scene.remove(planeScene);
      // scene.add(sphereScene);
    }
  }); 
  document.body.appendChild(sceneChangeButton);

  window.onresize = onWindowResize;
}

function initGui() {
  var gui = new dat.GUI();
  var dimensions = ["salinity", "temperature", "velocity"];
  gui
    .add({ Dimension: currentDimension }, "Dimension", dimensions)
    .onChange(function (value: string) {
      updateModelColor(value);
    });
}

function initRender() {
  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.setClearColor(0xffffff);
  renderer.outputEncoding = THREE.LinearEncoding;
  renderer.toneMapping = THREE.NoToneMapping;
  document.body.appendChild(renderer.domElement);
}

function initScene() {
  scene = new THREE.Scene();
}

function initCamera() {
  camera = new THREE.PerspectiveCamera(
    45,
    window.innerWidth / window.innerHeight,
    1,
    1000
  );
  camera.position.set(0, 0, 0);
  camera.lookAt(new THREE.Vector3(0, 0, 0));
}

function initLight() {
  scene.add(new THREE.AmbientLight(0x444444));
  light = new THREE.AmbientLight(0xffffff);
  light.position.set(0, 50, 50);
  light.castShadow = true;
  scene.add(light);
}

function initModel() {
  const helper = new THREE.AxesHelper(50);
  scene.add(helper);
  const loader = new THREE.FileLoader();
  loader.load('./ocean_data.json', function (jsonStr: any) {
    modelData = JSON.parse(jsonStr);
    //init Data
    data.salinity = modelData.salinity;
    data.temperature = modelData.temperature;
    data.velocity = vectorLength(modelData.velocity);
    ({ min: minSalinity, max: maxSalinity } = calculateDataRange(data.salinity, data.temperature, data.velocity, 'salinity'));
    ({ min: minTemperature, max: maxTemperature } = calculateDataRange(data.salinity, data.temperature, data.velocity, 'temperature'));
    ({ min: minVelocity, max: maxVelocity } = calculateDataRange(data.salinity, data.temperature, data.velocity, 'velocity'));
    createModelGeometry(modelData);
  })
}

function createModelGeometry(modelData: any) {
  const points: number[][] = modelData.points;
  const positions: Float32Array = new Float32Array(points.length * 3);
  const colors: Float32Array = new Float32Array(points.length * 3);
  //calculate the center point 
  const centerPoint: number[] = calculateCenterPoint(points);
  //calculate the positions
  for (let i = 0, j = 0; i < points.length; i++, j += 3) {
    //positions
    const point: number[] = points[i];
    positions[j] = point[0] - centerPoint[0];
    positions[j + 1] = point[1] - centerPoint[1];
    positions[j + 2] = point[2] - centerPoint[2];
  }
  //calculate the colors
  for (let i = 0, j = 0; i < points.length; i++, j += 3) {
    const currentValue: number = data.salinity[i];
    const colorValue: THREE.Color = getColorValue(currentValue, minSalinity, maxSalinity);
    colors[j] = colorValue.r;
    colors[j + 1] = colorValue.g;
    colors[j + 2] = colorValue.b;
  }
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
  const material = new THREE.PointsMaterial({ size: 0.3, vertexColors: true });
  const mesh: THREE.Points = new THREE.Points(geometry, material);
  scene.add(mesh);

  currentPoints = mesh;

  const boundingBox: THREE.Box3 = new THREE.Box3();
  boundingBox.setFromObject(mesh);
  const maxAxisSize: number = Math.max(boundingBox.max.x - boundingBox.min.x, boundingBox.max.y - boundingBox.min.y, boundingBox.max.z - boundingBox.min.z);
  const maxDistance: number = maxAxisSize / Math.tan((camera.fov / 2) * (Math.PI / 180));
  camera.position.set(0, 0, maxDistance * 1.5);
  camera.updateProjectionMatrix();
}

function vectorLength(vrctor: number[][]) {
  const length: number[] = [];
  for (let i = 0; i < vrctor.length; i++) {
    length.push(Math.sqrt(vrctor[i][0] * vrctor[i][0] + vrctor[i][1] * vrctor[i][1] + vrctor[i][2] * vrctor[i][2]));
  }
  return length;
}
function calculateCenterPoint(points: number[][]) {
  const centerPoint: number[] = [0, 0, 0];
  for (let i = 0; i < points.length; i++) {
    centerPoint[0] += points[i][0];
    centerPoint[1] += points[i][1];
    centerPoint[2] += points[i][2];
  }
  centerPoint[0] /= points.length;
  centerPoint[1] /= points.length;
  centerPoint[2] /= points.length;
  return centerPoint;
}

function calculateDataRange(salinity: number[], temperature: number[], velocity: number[], currentDimension: string) {
  let min: number = Infinity;
  let max: number = -Infinity;
  let value: number[];
  switch (currentDimension) {
    case "salinity":
      value = salinity;
      break;
    case "temperature":
      value = temperature;
      break;
    case "velocity":
      value = velocity;
      break;
    default:
      value = salinity;
  }
  for (let i = 0; i < value.length; i++) {
    if (value[i] < min) {
      min = value[i];
    }
    if (value[i] > max) {
      max = value[i];
    }
  }

  return { min: min, max: max };
}

function getColorValue(currentValue: number, minCurrentValue: number, maxCurrentValue: number) {
  const colorValue: THREE.Color = new THREE.Color();
  let normalizedValue: number = (currentValue - minCurrentValue) / (maxCurrentValue - minCurrentValue);
  if (normalizedValue < 0.5) {
    normalizedValue = 1 - normalizedValue * 2;
    colorValue.setRGB(
      colorMap.mi[0] * (1 - normalizedValue) + colorMap.lo[0] * normalizedValue,
      colorMap.mi[1] * (1 - normalizedValue) + colorMap.lo[1] * normalizedValue,
      colorMap.mi[2] * (1 - normalizedValue) + colorMap.lo[2] * normalizedValue);
  } else {
    normalizedValue = 1 - (normalizedValue - 0.5) * 2;
    colorValue.setRGB(
      colorMap.hi[0] * (1 - normalizedValue) + colorMap.mi[0] * normalizedValue,
      colorMap.hi[1] * (1 - normalizedValue) + colorMap.mi[1] * normalizedValue,
      colorMap.hi[2] * (1 - normalizedValue) + colorMap.mi[2] * normalizedValue);
  }
  return colorValue;
}
function initControls() {
  controls = new OrbitControls(camera, renderer.domElement);
  //controls.addEventListener( 'change', render );
  controls.enableDamping = true;
  controls.dampingFactor = 0.25;
  controls.enableZoom = true;
  controls.autoRotate = false;
  controls.autoRotateSpeed = 0.5;
  controls.minDistance = 1;
  controls.maxDistance = 200;
  controls.enablePan = true;
}

function initStats() {
  stats = new Stats();
  document.body.appendChild(stats.dom);
}

function animate() {
  render();
  stats.update();
  controls.update();
  requestAnimationFrame(animate);
}

function updateModelColor(dimension: string) {
  currentDimension = dimension;
  let updateData: number[] = [];
  let minCurrentValue: number = 0;
  let maxCurrentValue: number = 0;
  if (dimension === "velocity") {
    updateData = data.velocity;
    minCurrentValue = minVelocity;
    maxCurrentValue = maxVelocity;
  } else if (dimension === "salinity") {
    updateData = data.salinity;
    minCurrentValue = minSalinity;
    maxCurrentValue = maxSalinity;
  } else if (dimension === "temperature") {
    updateData = data.temperature;
    minCurrentValue = minTemperature;
    maxCurrentValue = maxTemperature;
  }

  const colors: Float32Array = new Float32Array(modelData.points.length * 3);
  for (let i = 0, j = 0; i < modelData.points.length; i++, j += 3) {
    const currentValue: number = updateData[i];
    const colorValue: THREE.Color = getColorValue(currentValue, minCurrentValue, maxCurrentValue);
    colors[j] = colorValue.r;
    colors[j + 1] = colorValue.g;
    colors[j + 2] = colorValue.b;
  }
  geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
  const material = new THREE.PointsMaterial({ size: 0.3, vertexColors: true });
  const mesh: THREE.Points = new THREE.Points(geometry, material);
  scene.add(mesh);

  if (currentPoints) {
    scene.remove(currentPoints);
    currentPoints.geometry.dispose();
  }
  currentPoints = mesh;
  render();
}

function render() {
  renderer.render(scene, camera);
}
// Function to handle scene change
function onSceneChange(sceneType: string) {
  if (sceneType === 'sphere') {
    currentScene = sphereScene;
    scene = sphereScene;
  } else if (sceneType === 'plane') {
    currentScene = planeScene;
    scene = planeScene;
  }
  updateModelColor(currentDimension);
  render();
}
function onWindowResize() {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  render();
  renderer.setSize(window.innerWidth, window.innerHeight);
}
