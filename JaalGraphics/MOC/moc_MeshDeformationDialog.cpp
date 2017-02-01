/****************************************************************************
** Meta object code from reading C++ file 'MeshDeformationDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshDeformationDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshDeformationDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshDeformationDialog_t {
    QByteArrayData data[15];
    char stringdata0[227];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshDeformationDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshDeformationDialog_t qt_meta_stringdata_JMeshDeformationDialog = {
    {
        QT_MOC_LITERAL(0, 0, 22), // "JMeshDeformationDialog"
        QT_MOC_LITERAL(1, 23, 9), // "resetData"
        QT_MOC_LITERAL(2, 33, 0), // ""
        QT_MOC_LITERAL(3, 34, 12), // "checkDisplay"
        QT_MOC_LITERAL(4, 47, 12), // "getImageMesh"
        QT_MOC_LITERAL(5, 60, 13), // "openLIMDialog"
        QT_MOC_LITERAL(6, 74, 25), // "openCurveShorteningDialog"
        QT_MOC_LITERAL(7, 100, 25), // "openMeshGeomQualityDialog"
        QT_MOC_LITERAL(8, 126, 17), // "getNewConstraints"
        QT_MOC_LITERAL(9, 144, 10), // "loadSource"
        QT_MOC_LITERAL(10, 155, 10), // "loadTarget"
        QT_MOC_LITERAL(11, 166, 9), // "runSolver"
        QT_MOC_LITERAL(12, 176, 21), // "setHandleAsConstraint"
        QT_MOC_LITERAL(13, 198, 16), // "setFixedBoundary"
        QT_MOC_LITERAL(14, 215, 11) // "closeDialog"

    },
    "JMeshDeformationDialog\0resetData\0\0"
    "checkDisplay\0getImageMesh\0openLIMDialog\0"
    "openCurveShorteningDialog\0"
    "openMeshGeomQualityDialog\0getNewConstraints\0"
    "loadSource\0loadTarget\0runSolver\0"
    "setHandleAsConstraint\0setFixedBoundary\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshDeformationDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    13,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   79,    2, 0x08 /* Private */,
    3,    0,   80,    2, 0x08 /* Private */,
    4,    0,   81,    2, 0x08 /* Private */,
    5,    0,   82,    2, 0x08 /* Private */,
    6,    0,   83,    2, 0x08 /* Private */,
    7,    0,   84,    2, 0x08 /* Private */,
    8,    0,   85,    2, 0x08 /* Private */,
    9,    0,   86,    2, 0x08 /* Private */,
    10,    0,   87,    2, 0x08 /* Private */,
    11,    0,   88,    2, 0x08 /* Private */,
    12,    0,   89,    2, 0x08 /* Private */,
    13,    0,   90,    2, 0x08 /* Private */,
    14,    0,   91,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMeshDeformationDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshDeformationDialog *_t = static_cast<JMeshDeformationDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->resetData();
            break;
        case 1:
            _t->checkDisplay();
            break;
        case 2:
            _t->getImageMesh();
            break;
        case 3:
            _t->openLIMDialog();
            break;
        case 4:
            _t->openCurveShorteningDialog();
            break;
        case 5:
            _t->openMeshGeomQualityDialog();
            break;
        case 6:
            _t->getNewConstraints();
            break;
        case 7:
            _t->loadSource();
            break;
        case 8:
            _t->loadTarget();
            break;
        case 9:
            _t->runSolver();
            break;
        case 10:
            _t->setHandleAsConstraint();
            break;
        case 11:
            _t->setFixedBoundary();
            break;
        case 12:
            _t->closeDialog();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshDeformationDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshDeformationDialog.data,
        qt_meta_data_JMeshDeformationDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshDeformationDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshDeformationDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshDeformationDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshDeformationDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshDeformationDialog"))
        return static_cast< Ui::MeshDeformationDialog*>(const_cast< JMeshDeformationDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshDeformationDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
